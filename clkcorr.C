#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2026 M. Keith, University of Manchester

/*
 *    This file is part of TEMPO2. 
 * 
 *    TEMPO2 is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    TEMPO2 is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing 
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */


 /** This is a re-write of the tempo2 clock correction logic that avoids the use of custom 
  * dynamic arrays and re-uses standard C++ infrastructure. It also avoids reading all files at startup
  * instead indexing the files and loading them as required 
  * 
  * Clock files are located in *.clk in the following paths:
  * * $TEMPO2/clock
  * * global variable tempo2_clock_path
  * 
  * Clock files files have a one-line header
  * # FROM TO [BADNESS]
  * where FROM is the source clock, e.g. UTC(JB)
  * TO is the target clock, e.g.  UTC(GPS)
  * Badness is an optional float that adds a penalty to the search path.
  * 
  * Then files contain lines with
  * MJD CORRECTION
  * where correction is in seconds. Lines starting with # should be ignored.
  * 
  * There may be multiple paths from clockFrom to clockTo and clock files don't always cover the full time range. E.g. there may be one route used for early data and one for later data.
  * When going from A to B, first build all routes from A to B and load the clock file contents. Then determine the time validity of that route, caching these so that future requests for that
  * observatory can re-use the same route. If no route is found, use Dijkstra's shortest path algorithm to find a route and load the files.
  * 
  * These functions are all exposed as C in the tempo2.h but are implemented in C++ and use the C++ standard library.
  * 
  * */


#include "tempo2.h"
#include "tabulatedfunction.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <glob.h>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct ClockCorrectionFunction {
  TabulatedFunction table;
  std::string fullPath;
  std::string fileName;
  std::string clockFrom;
  std::string clockTo;
  float badness;
  bool loaded;
  double startMJD;
  double endMJD;

  ClockCorrectionFunction()
    : badness(1.0f), loaded(false), startMJD(0.0), endMJD(0.0) {}
};

struct ClockSequence {
  std::vector<size_t> functionIndices;
  std::string clockFrom;
  std::string clockTo;
  std::string clockFromKey;
  std::string clockToKey;
  double startMJD;
  double endMJD;

  ClockSequence() : startMJD(0.0), endMJD(0.0) {}
};

struct PathBuildResult {
  std::vector<size_t> functionIndices;
  bool sourcePresent;
  bool targetPresent;
  bool reverseHint;

  PathBuildResult() : sourcePresent(false), targetPresent(false), reverseHint(false) {}
};

static bool g_clockCorrectionsInitialized = false;
static std::vector<ClockCorrectionFunction> g_clockCorrectionFunctions;
static std::vector<ClockSequence> g_clockCorrectionSequences;

static std::string trim(const std::string &s)
{
  const std::string ws(" \t\r\n");
  const std::string::size_type first = s.find_first_not_of(ws);
  if (first == std::string::npos)
    return std::string();
  const std::string::size_type last = s.find_last_not_of(ws);
  return s.substr(first, last - first + 1);
}

static std::string normalizeClockName(const std::string &name)
{
  std::string out = trim(name);
  std::transform(out.begin(), out.end(), out.begin(),
      [](unsigned char c){ return (char)std::toupper(c); });
  return out;
}

static bool clockNameEquals(const std::string &a, const std::string &b)
{
  return normalizeClockName(a) == normalizeClockName(b);
}

static std::string basenameFromPath(const std::string &path)
{
  const std::string::size_type pos = path.find_last_of('/');
  if (pos == std::string::npos)
    return path;
  return path.substr(pos + 1);
}

static bool parseClockHeaderLine(const std::string &header,
    std::string *clockFrom, std::string *clockTo, float *badness)
{
  char from[64];
  char to[64];
  float localBadness = 1.0f;

  if (header.empty() || header[0] != '#')
    return false;

  const int narg = std::sscanf(header.c_str() + 1, "%63s %63s %f", from, to, &localBadness);
  if (narg < 2)
    return false;

  *clockFrom = from;
  *clockTo = to;
  *badness = (narg == 3) ? localBadness : 1.0f;
  return true;
}

static void addIndexedClockFile(const std::string &filePath, int defaultBadnessOffset)
{
  FILE *f = std::fopen(filePath.c_str(), "r");
  if (f == NULL)
  {
    std::fprintf(stderr, "Fatal Error: Unable to open file %s for reading: %s\n",
        filePath.c_str(), std::strerror(errno));
    std::exit(1);
  }

  char headerLine[1024];
  if (std::fgets(headerLine, sizeof(headerLine), f) == NULL)
  {
    std::fprintf(stderr,
        "Error parsing clock file %s: file appears empty or corrupted.\n",
        filePath.c_str());
    std::fclose(f);
    std::exit(1);
  }
  std::fclose(f);

  std::string from;
  std::string to;
  float badness = 1.0f;
  if (!parseClockHeaderLine(headerLine, &from, &to, &badness))
  {
    std::fprintf(stderr,
        "Error parsing clock file %s: first line must be of form # clock_from clock_to [badness]\n",
        filePath.c_str());
    std::exit(1);
  }

  ClockCorrectionFunction func;
  func.fullPath = filePath;
  func.fileName = basenameFromPath(filePath);
  func.clockFrom = from;
  func.clockTo = to;
  func.badness = badness + (float)defaultBadnessOffset;
  g_clockCorrectionFunctions.push_back(func);
}

static std::vector<std::string> globClockFiles(const std::string &pattern)
{
  glob_t g;
  std::vector<std::string> files;

  const int globRet = glob(pattern.c_str(), 0, NULL, &g);
  if (globRet == GLOB_NOSPACE)
  {
    std::fprintf(stderr, "Out of memory in clkcorr.C\n");
    std::exit(1);
  }
#ifdef GLOB_ABORTED
  if (globRet == GLOB_ABORTED)
  {
    std::fprintf(stderr, "Read error in clkcorr.C\n");
    std::exit(1);
  }
#endif
  if (globRet == 0)
  {
    for (char **pfname = g.gl_pathv; *pfname != NULL; ++pfname)
      files.push_back(*pfname);
  }
  globfree(&g);

  std::sort(files.begin(), files.end());
  return files;
}

static void ensureClockCorrectionsInitialized(int dispWarnings)
{
  if (g_clockCorrectionsInitialized)
    return;

  int defaultBadnessOffset = 0;

  if (std::strlen(tempo2_clock_path) > 0)
  {
    logmsg("Note, using '%s' to look for extra clock files", tempo2_clock_path);
    const std::string pattern = std::string(tempo2_clock_path) + "/*.clk";
    const std::vector<std::string> extraFiles = globClockFiles(pattern);

    if (extraFiles.empty())
      logerr("No clock correction files in '%s'", tempo2_clock_path);
    else
    {
      for (const std::string &path : extraFiles)
      {
        addIndexedClockFile(path, 0);
        const ClockCorrectionFunction &func = g_clockCorrectionFunctions.back();
        logmsg("Loaded extra clock correction %s : %s -> %s badness %.3g",
            func.fileName.c_str(), func.clockFrom.c_str(), func.clockTo.c_str(), func.badness);
      }
      defaultBadnessOffset = 1;
    }
  }

  const char *tempo2 = std::getenv(TEMPO2_ENVIRON);
  if (tempo2 == NULL || std::strlen(tempo2) == 0)
  {
    displayMsg(2, "CLK11", "TEMPO2 environment path is not set", "", dispWarnings);
    std::exit(1);
  }

  const std::string defaultPattern = std::string(tempo2) + "/clock/*.clk";
  const std::vector<std::string> defaultFiles = globClockFiles(defaultPattern);
  if (defaultFiles.empty())
  {
    std::fprintf(stderr, "No clock correction files in $TEMPO2/clock\n");
    std::exit(1);
  }

  for (const std::string &path : defaultFiles)
  {
    addIndexedClockFile(path, defaultBadnessOffset);
    const ClockCorrectionFunction &func = g_clockCorrectionFunctions.back();
    logdbg("Loaded clock correction %s : %s -> %s badness %.3g",
        func.fileName.c_str(), func.clockFrom.c_str(), func.clockTo.c_str(), func.badness);
  }

  g_clockCorrectionsInitialized = true;
}

static void ensureFunctionLoaded(ClockCorrectionFunction *func)
{
  if (func->loaded)
    return;

  TabulatedFunction_load(&func->table, const_cast<char *>(func->fullPath.c_str()));
  func->startMJD = TabulatedFunction_getStartX(&func->table);
  func->endMJD = TabulatedFunction_getEndX(&func->table);
  func->loaded = true;
}

static bool functionValidAtMJD(ClockCorrectionFunction *func, double mjd)
{
  ensureFunctionLoaded(func);
  return func->startMJD <= mjd && func->endMJD >= mjd;
}

static std::pair<double,double> sequenceMJDRange(const std::vector<size_t> &indices)
{
  double start = 0.0;
  double end = 1e6;

  for (size_t idx : indices)
  {
    ClockCorrectionFunction *func = &g_clockCorrectionFunctions[idx];
    ensureFunctionLoaded(func);
    if (func->startMJD > start)
      start = func->startMJD;
    if (func->endMJD < end)
      end = func->endMJD;
  }

  return std::make_pair(start, end);
}

static bool sequenceMatches(const ClockSequence &seq,
    const std::string &clockFrom, const std::string &clockTo, double mjd)
{
  return seq.clockFromKey == normalizeClockName(clockFrom)
    && seq.clockToKey == normalizeClockName(clockTo)
    && seq.startMJD <= mjd
    && seq.endMJD >= mjd;
}

static int findFunctionIndexFromToken(const std::string &token)
{
  for (size_t i = 0; i < g_clockCorrectionFunctions.size(); ++i)
  {
    if (g_clockCorrectionFunctions[i].fileName == token)
      return (int)i;
  }

  const std::string tokenBase = basenameFromPath(token);
  for (size_t i = 0; i < g_clockCorrectionFunctions.size(); ++i)
  {
    if (g_clockCorrectionFunctions[i].fileName == tokenBase)
      return (int)i;
  }

  return -1;
}

static PathBuildResult buildDirectedPath(const std::string &clockFrom,
    const std::string &clockTo, double mjd)
{
  PathBuildResult result;

  struct Edge {
    int from;
    int to;
    size_t funcIndex;
    float weight;
  };

  std::vector<std::string> nodeKeys;
  std::unordered_map<std::string, int> nodeIndex;
  std::vector<Edge> edges;

  auto getNode = [&nodeKeys, &nodeIndex](const std::string &name) -> int {
    const std::string key = normalizeClockName(name);
    std::unordered_map<std::string, int>::iterator it = nodeIndex.find(key);
    if (it != nodeIndex.end())
      return it->second;

    const int idx = (int)nodeKeys.size();
    nodeKeys.push_back(key);
    nodeIndex[key] = idx;
    return idx;
  };

  const std::string fromKey = normalizeClockName(clockFrom);
  const std::string toKey = normalizeClockName(clockTo);

  for (size_t i = 0; i < g_clockCorrectionFunctions.size(); ++i)
  {
    ClockCorrectionFunction *func = &g_clockCorrectionFunctions[i];
    if (!functionValidAtMJD(func, mjd))
      continue;

    const int fromNode = getNode(func->clockFrom);
    const int toNode = getNode(func->clockTo);
    Edge e;
    e.from = fromNode;
    e.to = toNode;
    e.funcIndex = i;
    e.weight = func->badness;
    edges.push_back(e);
  }

  result.sourcePresent = (nodeIndex.find(fromKey) != nodeIndex.end());
  result.targetPresent = (nodeIndex.find(toKey) != nodeIndex.end());

  if (!result.sourcePresent || !result.targetPresent)
    return result;

  const int source = nodeIndex[fromKey];
  const int target = nodeIndex[toKey];

  std::vector<std::vector<int> > adjacency(nodeKeys.size());
  std::vector<std::vector<int> > reverseAdjacency(nodeKeys.size());
  for (size_t i = 0; i < edges.size(); ++i)
  {
    adjacency[edges[i].from].push_back((int)i);
    reverseAdjacency[edges[i].to].push_back(edges[i].from);
  }

  std::vector<float> dist(nodeKeys.size(), FLT_MAX);
  std::vector<int> visited(nodeKeys.size(), 0);
  std::vector<int> prevEdge(nodeKeys.size(), -1);
  std::vector<int> prevNode(nodeKeys.size(), -1);
  dist[source] = 0.0f;

  while (true)
  {
    int bestNode = -1;
    float best = FLT_MAX;
    for (size_t i = 0; i < dist.size(); ++i)
    {
      if (!visited[i] && dist[i] < best)
      {
        best = dist[i];
        bestNode = (int)i;
      }
    }

    if (bestNode < 0)
      break;
    visited[bestNode] = 1;
    if (bestNode == target)
      break;

    for (int edgeIdx : adjacency[bestNode])
    {
      const Edge &edge = edges[edgeIdx];
      const float newDist = dist[bestNode] + edge.weight;
      if (newDist < dist[edge.to])
      {
        dist[edge.to] = newDist;
        prevEdge[edge.to] = edgeIdx;
        prevNode[edge.to] = bestNode;
      }
    }
  }

  if (!visited[target])
  {
    if (source < (int)reverseAdjacency.size() && target < (int)reverseAdjacency.size())
    {
      std::vector<int> stack;
      std::vector<int> seen(nodeKeys.size(), 0);
      stack.push_back(source);
      seen[source] = 1;
      while (!stack.empty())
      {
        const int v = stack.back();
        stack.pop_back();
        if (v == target)
        {
          result.reverseHint = true;
          break;
        }
        for (int next : reverseAdjacency[v])
        {
          if (!seen[next])
          {
            seen[next] = 1;
            stack.push_back(next);
          }
        }
      }
    }
    return result;
  }

  std::vector<size_t> reversed;
  for (int v = target; v != source; v = prevNode[v])
  {
    if (v < 0 || prevEdge[v] < 0)
    {
      reversed.clear();
      break;
    }
    reversed.push_back(edges[prevEdge[v]].funcIndex);
  }
  std::reverse(reversed.begin(), reversed.end());
  result.functionIndices.swap(reversed);
  return result;
}

static const ClockSequence *cacheSequence(const std::vector<size_t> &indices,
    const std::string &clockFrom, const std::string &clockTo)
{
  ClockSequence seq;
  seq.functionIndices = indices;
  seq.clockFrom = clockFrom;
  seq.clockTo = clockTo;
  seq.clockFromKey = normalizeClockName(clockFrom);
  seq.clockToKey = normalizeClockName(clockTo);

  const std::pair<double,double> range = sequenceMJDRange(indices);
  seq.startMJD = range.first;
  seq.endMJD = range.second;

  g_clockCorrectionSequences.push_back(seq);
  return &g_clockCorrectionSequences.back();
}

static const ClockSequence *getClockCorrectionSequenceInternal(const std::string &clockFrom,
    const std::string &clockTo, double mjd, int warnings)
{
  ensureClockCorrectionsInitialized(warnings);

  for (const ClockSequence &seq : g_clockCorrectionSequences)
  {
    if (sequenceMatches(seq, clockFrom, clockTo, mjd))
      return &seq;
  }

  logdbg("Making clock sequence from %s to %s", clockFrom.c_str(), clockTo.c_str());

  PathBuildResult result = buildDirectedPath(clockFrom, clockTo, mjd);
  if (!result.sourcePresent)
  {
    char msg[1000], msg2[1000];
    std::sprintf(msg, "no clock corrections available for clock %s for MJD", clockFrom.c_str());
    std::sprintf(msg2, "%.1f", mjd);
    displayMsg(2, "CLK3", msg, msg2, warnings);
    return NULL;
  }
  if (!result.targetPresent)
  {
    char msg[1000], msg2[1000];
    std::sprintf(msg, "no clock corrections available for clock %s for MJD", clockTo.c_str());
    std::sprintf(msg2, "%.1f", mjd);
    displayMsg(2, "CLK3", msg, msg2, warnings);
    return NULL;
  }
  if (result.functionIndices.empty())
  {
    char msg[1000], msg2[1000];
    std::sprintf(msg, "no directed clock corrections available from %s to %s for MJD",
        clockFrom.c_str(), clockTo.c_str());
    if (result.reverseHint)
      std::sprintf(msg2, "%.1f (reverse path exists but reverse traversal is disabled)", mjd);
    else
      std::sprintf(msg2, "%.1f", mjd);
    displayMsg(2, "CLK7", msg, msg2, warnings);
    return NULL;
  }

  if (warnings == 0)
  {
    std::printf("Using the following chain of clock corrections for %s -> %s\n",
        clockFrom.c_str(), clockTo.c_str());
    for (size_t idx : result.functionIndices)
    {
      const ClockCorrectionFunction &func = g_clockCorrectionFunctions[idx];
      std::printf("%s : %s -> %s\n",
          func.fileName.c_str(), func.clockFrom.c_str(), func.clockTo.c_str());
    }
  }

  return cacheSequence(result.functionIndices, clockFrom, clockTo);
}

static bool fillCorrectionsFromSequence(observation *obs,
    const ClockSequence *sequence, const std::string &clockTo,
    std::string currentClock, double *totalCorrection)
{
  for (size_t ifunc = 0;
      ifunc < sequence->functionIndices.size() && !clockNameEquals(currentClock, clockTo);
      ++ifunc)
  {
    if (obs->nclock_correction >= MAX_CLK_CORR)
    {
      displayMsg(2, "CLK12", "clock correction chain exceeded MAX_CLK_CORR", "", 0);
      return false;
    }

    ClockCorrectionFunction *func = &g_clockCorrectionFunctions[sequence->functionIndices[ifunc]];
    ensureFunctionLoaded(func);

    if (!clockNameEquals(currentClock, func->clockFrom))
    {
      logerr("Broken directed clock correction chain: expected %s, got %s",
          currentClock.c_str(), func->clockFrom.c_str());
      return false;
    }

    const double stepCorrection = TabulatedFunction_getValue(&func->table,
        (double)obs->sat + *totalCorrection / SECDAY);

    obs->correctionsTT[obs->nclock_correction].correction = stepCorrection;
    std::snprintf(obs->correctionsTT[obs->nclock_correction].corrects_to,
        sizeof(obs->correctionsTT[obs->nclock_correction].corrects_to),
        "%s", func->clockTo.c_str());
    ++obs->nclock_correction;

    *totalCorrection += stepCorrection;
    currentClock = func->clockTo;
  }

  return clockNameEquals(currentClock, clockTo);
}

static bool accumulateCorrectionFromSequence(const ClockSequence *sequence,
    observation *obs, const std::string &clockTo,
    std::string currentClock, double *correction)
{
  for (size_t ifunc = 0;
      ifunc < sequence->functionIndices.size() && !clockNameEquals(currentClock, clockTo);
      ++ifunc)
  {
    ClockCorrectionFunction *func = &g_clockCorrectionFunctions[sequence->functionIndices[ifunc]];
    ensureFunctionLoaded(func);

    if (!clockNameEquals(currentClock, func->clockFrom))
    {
      logerr("Broken directed clock correction chain: expected %s, got %s",
          currentClock.c_str(), func->clockFrom.c_str());
      return false;
    }

    *correction += TabulatedFunction_getValue(&func->table,
        (double)obs->sat + *correction / SECDAY);
    currentClock = func->clockTo;
  }

  return clockNameEquals(currentClock, clockTo);
}

} // namespace




/* defineClockCorrectionSequence: call to provide the clock correction
module with a sequence of files to use for corrections. May be called
multiple times for sequences with different start/end clocks (e.g. for
multi-observatory fitting). */
void defineClockCorrectionSequence(char *fileList,int dispWarnings){
  ensureClockCorrectionsInitialized(dispWarnings);

  std::stringstream ss(fileList == NULL ? "" : fileList);
  std::string token;
  std::vector<size_t> indices;

  while (ss >> token)
  {
    const int idx = findFunctionIndexFromToken(token);
    if (idx < 0)
    {
      std::printf("Requested clock correction file %s not found!\n", token.c_str());
      std::exit(1);
    }
    indices.push_back((size_t)idx);
  }

  if (indices.empty())
  {
    displayMsg(2, "CLK13", "Empty CLK_CORR_CHAIN definition", "", dispWarnings);
    std::exit(1);
  }

  for (size_t i = 1; i < indices.size(); ++i)
  {
    const ClockCorrectionFunction &prev = g_clockCorrectionFunctions[indices[i - 1]];
    const ClockCorrectionFunction &next = g_clockCorrectionFunctions[indices[i]];
    if (!clockNameEquals(prev.clockTo, next.clockFrom))
    {
      std::fprintf(stderr,
          "Invalid directed clock correction chain: %s (%s -> %s) cannot precede %s (%s -> %s)\n",
          prev.fileName.c_str(), prev.clockFrom.c_str(), prev.clockTo.c_str(),
          next.fileName.c_str(), next.clockFrom.c_str(), next.clockTo.c_str());
      std::exit(1);
    }
  }

  const ClockCorrectionFunction &first = g_clockCorrectionFunctions[indices.front()];
  const ClockCorrectionFunction &last = g_clockCorrectionFunctions[indices.back()];
  cacheSequence(indices, first.clockFrom, last.clockTo);
}

/* getClockCorrections : gets the sequence of corrections for a particular
observation and stores them in obs->clock_corrections.  Uses one
of the pre-defined sequences (e.g. from defineClockCorrectionSequence) if
available, otherwise makes one automatically.  */
void getClockCorrections(observation *obs, const char *clockFrom, const char *clockTo, int warnings){
  obs->nclock_correction = 0;

  std::string fromClock = (clockFrom == NULL ? "" : clockFrom);
  std::string targetClock = (clockTo == NULL ? "" : clockTo);
  if (fromClock.empty())
    fromClock = getObservatory(obs->telID)->clock_name;

  const std::string originalFrom = fromClock;
  if (clockNameEquals(fromClock, targetClock))
    return;

  const ClockSequence *sequence = getClockCorrectionSequenceInternal(fromClock, targetClock,
      (double)obs->sat, warnings);
  bool usedApproximation = false;

  if (sequence == NULL)
  {
    char msg[1000], msg2[1000];
    std::sprintf(msg, "Trying assuming UTC =");
    std::sprintf(msg2, "%s", fromClock.c_str());
    displayMsg(1, "CLK4", msg, msg2, warnings);

    fromClock = "UTC";
    if (clockNameEquals(fromClock, targetClock))
      return;

    sequence = getClockCorrectionSequenceInternal(fromClock, targetClock,
        (double)obs->sat, warnings);
    usedApproximation = true;
    if (sequence == NULL)
    {
      std::sprintf(msg, "Trying TT(TAI) instead of ");
      std::sprintf(msg2, "%s", targetClock.c_str());
      displayMsg(2, "CLK5", msg, msg2, warnings);

      fromClock = originalFrom;
      targetClock = "TT(TAI)";
      if (clockNameEquals(fromClock, targetClock))
        return;

      sequence = getClockCorrectionSequenceInternal(fromClock, targetClock,
          (double)obs->sat, warnings);
      if (sequence == NULL)
      {
        displayMsg(2, "CLK8", "Trying both", "", warnings);

        fromClock = "UTC";
        if (clockNameEquals(fromClock, targetClock))
          return;

        sequence = getClockCorrectionSequenceInternal(fromClock, targetClock,
            (double)obs->sat, warnings);
        if (sequence == NULL)
        {
          if (warnings == 0)
          {
            std::printf("Warning [CLK:7], no directed clock correction available for TOA @ MJD %.4lf!\n",
                (double)obs->sat);
          }
          return;
        }
      }
    }
  }

  double correction = 0.0;
  if (!fillCorrectionsFromSequence(obs, sequence, targetClock, fromClock, &correction))
  {
    displayMsg(2, "CLK14", "Broken directed clock correction chain", "", warnings);
    obs->nclock_correction = 0;
    return;
  }

  if (usedApproximation)
  {
    displayMsg(1, "CLK9", "... ok, using stated approximation", "", warnings);
  }
}

/* getCorrectionTT : convenience function to return the sum of all
correctionsTT terms in an observation */
double getCorrectionTT(observation *obs){
  double correction = 0.0;
  for (int ic = 0; ic < obs->nclock_correction; ++ic)
    correction += obs->correctionsTT[ic].correction;
  return correction;
}
/* convenience function to obtain correction to a named clock 
(for intermediate use e.g. in obtaining Earth orientation parameters;
does not store steps used in obs->correctionsTT */
double getCorrection(observation *obs, const char *clockFrom, const char *clockTo, int warnings){
  observatory *site = NULL;
  const char *CVS_verNum = "$Id$";

  std::string fromClock = (clockFrom == NULL ? "" : clockFrom);
  std::string targetClock = (clockTo == NULL ? "" : clockTo);
  if (fromClock.empty())
    site = getObservatory(obs->telID);

  if (displayCVSversion == 1)
    CVSdisplayVersion("clkcorr.C", "getCorrection()", CVS_verNum);

  if (fromClock.empty() && site != NULL)
    fromClock = site->clock_name;

  if (clockNameEquals(fromClock, targetClock))
    return 0.0;

  std::string currentClock = fromClock;
  const ClockSequence *sequence = getClockCorrectionSequenceInternal(fromClock,
      targetClock, (double)obs->sat, warnings);

  if (sequence == NULL)
  {
    char msg[1000], msg2[1000];
    std::sprintf(msg, "Proceeding assuming UTC = ");
    std::sprintf(msg2, "%s", currentClock.c_str());
    displayMsg(1, "CLK6", msg, msg2, warnings);

    if (clockNameEquals(targetClock, "UTC"))
      return 0.0;

    sequence = getClockCorrectionSequenceInternal("UTC", targetClock,
        (double)obs->sat, warnings);
    currentClock = "UTC";
    if (sequence == NULL)
    {
      if (warnings == 0)
        std::printf("!Warning [CLK:9], no directed clock correction available for TOA @ MJD %.4lf!\n",
            (double)obs->sat);
      return 0.0;
    }
  }

  double correction = 0.0;
  if (!accumulateCorrectionFromSequence(sequence, obs, targetClock, currentClock, &correction))
  {
    displayMsg(2, "CLK14", "Broken directed clock correction chain", "", warnings);
    return 0.0;
  }

  return correction;
}