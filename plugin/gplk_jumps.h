/** JUMP display names and fixed instrumental corrections. GPLv3-or-later.
 * This policy belongs to gplk; it does not change TEMPO2's parameter parser or
 * edit the input par file. Fixing a correction clears only its fit flag, never
 * its value, uncertainty, SAT convention, selector or observation membership.
 */
#ifndef TEMPO2_GPLK_JUMPS_H
#define TEMPO2_GPLK_JUMPS_H
#include "tempo2.h"
#include <string>
#include <sstream>
#include <cctype>

struct GplkJumpDescription {
    std::string selector,match,label;
    bool correction;
    GplkJumpDescription():correction(false) {}
};

/** Conservative convention: -<instrument>_<5-digit MJD>_jump 1.
 * Do not treat all numeric selectors, all non-sys flags, or all *_jump names
 * as corrections: those can be legitimate user-selectable fit parameters. */
static bool gplkIsCorrectionSelector(const std::string &selector,const std::string &match)
{
    if (match!="1" || selector.empty() || selector[0]!='-') return false;
    std::string name=selector;
    for (std::size_t i=0;i<name.size();++i)
        name[i]=static_cast<char>(std::tolower(static_cast<unsigned char>(name[i])));
    const std::string suffix="_jump";
    if (name.size()<suffix.size()+8 || name.compare(name.size()-suffix.size(),suffix.size(),suffix)!=0)
        return false;
    const std::size_t end=name.size()-suffix.size(),start=end-5;
    if (name[start-1]!='_' || start<3) return false;
    for (std::size_t i=start;i<end;++i)
        if (!std::isdigit(static_cast<unsigned char>(name[i]))) return false;
    return true;
}

static GplkJumpDescription gplkDescribeJump(const char *raw)
{
    GplkJumpDescription result;
    std::istringstream input(raw);
    input>>result.selector>>result.match;
    result.correction=gplkIsCorrectionSelector(result.selector,result.match);
    if (!result.selector.empty() && result.selector[0]=='-') {
        // Conventional grouping flags use their selected value as the name.
        // Dedicated boolean flags carry the jump's name in the selector itself:
        // "-blahblah 1 <offset>" must display "blahblah", never just "1".
        const bool group=result.selector=="-sys" || result.selector=="-f" ||
                         result.selector=="-g" || result.selector=="-be" ||
                         result.selector=="-backend" || result.selector=="-group";
        const bool namedFlag=result.match=="1" && !group;
        result.label=result.match.empty() || result.correction || namedFlag ?
                     result.selector.substr(1) : result.match;
    } else if (result.selector.empty()) {
        result.label="Unnamed jump";
    } else {
        result.label=result.selector+" "+result.match;
        if (result.selector=="MJD" || result.selector=="FREQ") {
            std::string end; if (input>>end) result.label+=".."+end;
        }
    }
    return result;
}

static int gplkKeepCorrectionsFixed(pulsar *psr,int count)
{
    int changed=0;
    for (int p=0;p<count;++p) for (int i=1;i<=psr[p].nJumps;++i)
        if (psr[p].fitJump[i]!=0 && gplkDescribeJump(psr[p].jumpStr[i]).correction) {
            psr[p].fitJump[i]=0; ++changed;
        }
    return changed;
}

static bool gplkToggleJumpFit(pulsar &p,int index)
{
    if (index<1 || index>p.nJumps || gplkDescribeJump(p.jumpStr[index]).correction) return false;
    p.fitJump[index]=!p.fitJump[index];
    if (!p.fitJump[index]) p.jumpValErr[index]=0;
    return true;
}
#endif
