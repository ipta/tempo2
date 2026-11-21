/**
 * @file
 * @brief Native Giza devices, appearance, overlays and input for gplk.
 *
 * Helpers borrow plot arrays; they never change observations or fit results.
 * Coordinates are in the caller's current world window unless stated otherwise.
 * Drawing helpers restore attributes, but Giza save/restore does not restore
 * viewport/window geometry, selected devices or colour-table definitions.
 *
 * Internal-linkage state belongs to the single plugin translation unit. See
 * documentation/GPLK_GIZA_DEVELOPER_GUIDE.md before splitting this into .C files.
 * Distributed under TEMPO2's GPLv3-or-later licence.
 */
#ifndef TEMPO2_GPLK_GIZA_H
#define TEMPO2_GPLK_GIZA_H

extern "C" {
#include <giza.h>
}
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

enum GplkTheme { GPLK_DARK, GPLK_LIGHT, GPLK_CLASSIC };
enum GplkPalette { GPLK_OCEAN, GPLK_CONTRAST, GPLK_SUNSET };
static GplkTheme gplkTheme = GPLK_DARK;
static GplkPalette gplkPalette = GPLK_OCEAN;
static bool gplkThemeExplicit = false, gplkPaletteExplicit = false;
static float gplkPointScale = 1.3f;
// All rectangles use full-surface, bottom-left normalised coordinates. Keep
// drawing and hit testing on these same rectangles when changing the layout.
struct GplkRect {
    float left,right,bottom,top;
    GplkRect(float l=0,float r=1,float b=0,float t=1):left(l),right(r),bottom(b),top(t) {}
    bool contains(float x,float y) const { return x>=left && x<right && y>=bottom && y<top; }
    float localX(float x) const { return (x-left)/(right-left); }
    float localY(float y) const { return (y-bottom)/(top-bottom); }
    void select() const { giza_set_viewport_float(left,right,bottom,top); }
};
struct GplkLayout {
    GplkRect top,side,bottom,jumps,colours,plot;
    float width,height,uiFont,sideFont,plotFont;
    int columns,visibleRows;
    void update(float w,float h,int entries) {
        width=std::max(w,1.0f); height=std::max(h,1.0f);
        columns=std::max(1,std::min(8,int(width/190)));
        const int rows=(entries+columns-1)/columns;
        const float topPixels=std::min(112.0f,height*0.20f);
        const float bottomPixels=std::min(82.0f+24*rows,height*0.42f);
        visibleRows=std::max(0,int((bottomPixels-82)/24));
        top=GplkRect(0,1,1-topPixels/height,1);
        bottom=GplkRect(0,1,0,bottomPixels/height);
        jumps=GplkRect(0,1,0,std::max(1.0f,bottomPixels-56)/height);
        const float sidebar=std::min(225.0f,width*0.25f);
        const float marginX=std::min(65.0f,width*0.075f);
        const float marginY=std::min(48.0f,height*0.075f);
        const float colourWidth=std::min(245.0f,width*0.23f);
        colours=GplkRect(1-colourWidth/width,1,bottom.top,top.bottom);
        plot=GplkRect((sidebar+marginX)/width,colours.left-marginX/width,
                     bottom.top+marginY/height,top.bottom-marginY/height);
        side=GplkRect(0,sidebar/width,plot.bottom,plot.top);
        plotFont=600/height;
        uiFont=plotFont;
        sideFont=plotFont*std::min(1.0f,std::min(sidebar/225,
                                      (plot.top-plot.bottom)*height/260));
    }
};
static GplkLayout gplkLayout;
static bool gplkResponsiveFonts=false;
static float gplkPlotCharacterHeight(float size) { return size*(gplkResponsiveFonts ? gplkLayout.plotFont : 1); }
static float gplkUiCharacterHeight(float size) { return size*(gplkResponsiveFonts ? gplkLayout.uiFont : 1); }

/** Use the double native batch entry point. Some Giza builds unconditionally
 * end the caller's frame buffer in giza_points_float; giza_points preserves it.
 * Reuse conversion storage across groups/frames instead of allocating per TOA. */
static void gplkDrawPoints(int n,const float *x,const float *y,int symbol)
{
    if (n<=0) return;
    static std::vector<double> dx,dy;
    dx.resize(n); dy.resize(n);
    for (int i=0;i<n;++i) { dx[i]=x[i]; dy[i]=y[i]; }
    giza_points(n,&dx[0],&dy[0],symbol);
}

static void gplkTextSize(const char *text,float *width,float *height)
{
    // Giza's WORLD text-length conversion is incorrect in some versions.
    // Query pixels and convert distances through the actual viewport/window.
    float vx1,vx2,vy1,vy2,wx1,wx2,wy1,wy2,pw,ph;
    giza_qtextlen_float(GIZA_UNITS_PIXELS,text,width,height);
    giza_get_viewport_float(GIZA_UNITS_NORMALIZED,&vx1,&vx2,&vy1,&vy2);
    giza_get_window_float(&wx1,&wx2,&wy1,&wy2);
    giza_get_paper_size_float(GIZA_UNITS_PIXELS,&pw,&ph);
    *width*=std::fabs((wx2-wx1)/((vx2-vx1)*pw));
    *height*=std::fabs((wy2-wy1)/((vy2-vy1)*ph));
}

/** UI identifiers are literal, whereas Giza text treats _ ^ { } and backslash
 * as mathematical markup. Draw ordinary runs together and the few special
 * glyphs separately so backend_01 is neither modified nor shown as a subscript.
 * Return width in current world coordinates; paint=false measures only. */
static float gplkLiteralText(float x,float y,const std::string &text,bool paint=true)
{
    const float start=x;
    for (std::size_t at=0;at<text.size();) {
        const std::size_t special=text.find_first_of("_^{}\\",at);
        const std::size_t end=special==std::string::npos ? text.size() : special;
        if (end>at) {
            const std::string run=text.substr(at,end-at);
            float w,h; gplkTextSize(run.c_str(),&w,&h);
            if (paint) giza_text_float(x,y,run.c_str());
            x+=w;
        }
        if (end==text.size()) break;
        const char c=text[end];
        float w,h; gplkTextSize("n",&w,&h);
        if (c=='{' || c=='}') {
            const char *glyph=c=='{' ? "\\(1225)" : "\\(1226)";
            gplkTextSize(glyph,&w,&h);
            if (paint) giza_text_float(x,y,glyph);
        } else if (paint) {
            giza_save(); giza_set_line_style(GIZA_LS_SOLID); giza_set_line_width(.8);
            if (c=='_') { giza_move_float(x,y-h*.12f); giza_draw_float(x+w,y-h*.12f); }
            else if (c=='^') {
                giza_move_float(x,y+h*.65f); giza_draw_float(x+w*.5f,y+h); giza_draw_float(x+w,y+h*.65f);
            } else { giza_move_float(x,y+h); giza_draw_float(x+w,y); }
            giza_restore();
        }
        x+=w; at=end+1;
    }
    return x-start;
}
// Reserve high indices for UI/scratch/overlay roles; 0..15 retain data meanings.
// A palette change must not turn a bin or model overlay into a residual series.
enum { GPLK_PANEL = 256, GPLK_GRID, GPLK_MUTED, GPLK_ACCENT, GPLK_ERROR,
       GPLK_BINS, GPLK_MODEL, GPLK_SELECTED, GPLK_SOLAR, GPLK_NOISE_OVERLAY };

/** Parse configuration only; apply it with gplkApplyTheme() on an open device. */
static bool gplkSetTheme(const char *name)
{
    if (std::strcmp(name, "dark") == 0) gplkTheme = GPLK_DARK;
    else if (std::strcmp(name, "light") == 0) gplkTheme = GPLK_LIGHT;
    else if (std::strcmp(name, "classic") == 0) gplkTheme = GPLK_CLASSIC;
    else return false;
    return true;
}

static bool gplkSetPalette(const char *name)
{
    if (std::strcmp(name, "ocean") == 0) gplkPalette = GPLK_OCEAN;
    else if (std::strcmp(name, "contrast") == 0) gplkPalette = GPLK_CONTRAST;
    else if (std::strcmp(name, "sunset") == 0) gplkPalette = GPLK_SUNSET;
    else return false;
    return true;
}

static void gplkSetRGB(int index, unsigned int rgb)
{
    giza_set_colour_representation(index, ((rgb >> 16) & 255) / 255.0,
                                   ((rgb >> 8) & 255) / 255.0, (rgb & 255) / 255.0);
}

/** Apply theme/series colours to the current device; preserves plot geometry. */
static void gplkApplyTheme()
{
    if (gplkTheme == GPLK_CLASSIC) return;
    const bool dark = gplkTheme == GPLK_DARK;
    gplkSetRGB(0, dark ? 0x111822 : 0xf3f6fa);
    gplkSetRGB(1, dark ? 0xdce6f2 : 0x263749);
    gplkSetRGB(GPLK_PANEL, dark ? 0x192331 : 0xffffff);
    gplkSetRGB(GPLK_GRID, dark ? 0x2c3b4e : 0xe1e8f0);
    gplkSetRGB(GPLK_MUTED, dark ? 0x91a3b7 : 0x66798d);
    // [GplkPalette enum][0=light, 1=dark][six series colours]. Keep enum/parser,
    // toolbar names and its cycle modulus in sync when adding a palette.
    static const unsigned int colours[3][2][6] = {
        {{0x2463a6,0x147d70,0x9c660f,0x7754ae,0xbe4667,0x247c92},
         {0x70b8ff,0x67d5b5,0xf2c66d,0xb7a0ef,0xf58d9b,0x6bd0d6}},
        {{0x0072b2,0xd55e00,0x007f5f,0xa64e88,0x936600,0x26728c},
         {0x69b9ec,0xf2a56f,0x68d3af,0xe1a4ce,0xf0d36c,0x87d5e8}},
        {{0x9a4168,0xb76622,0x6652a2,0x237f88,0xb13e45,0x537d32},
         {0xe99bb9,0xf4bb82,0xb8a5ef,0x7fd2d8,0xf08d94,0xb7d78c}}
    };
    for (int i = 2; i <= 12; ++i)
        gplkSetRGB(i, colours[gplkPalette][dark ? 1 : 0][(i-2)%6]);
    gplkSetRGB(13, dark ? 0xf58d9b : 0xbe4667);
    gplkSetRGB(14, dark ? 0x91a3b7 : 0x66798d);
    gplkSetRGB(15, dark ? 0x67788e : 0x9aabba);
    gplkSetRGB(GPLK_ACCENT, colours[gplkPalette][dark ? 1 : 0][0]);
}

/** Fill a rectangle in current world coordinates without changing the window. */
static void gplkDrawPanel(float x1, float x2, float y1, float y2)
{
    if (gplkTheme == GPLK_CLASSIC) return;
    giza_save();
    giza_set_colour_index(GPLK_PANEL);
    giza_set_fill(GIZA_FILL_SOLID);
    giza_rectangle_float(x1, x2, y1, y2);
    giza_restore();
}

static void gplkDrawPlotBackground(float x1, float x2, float y1, float y2)
{
    if (gplkTheme == GPLK_CLASSIC) return;
    gplkDrawPanel(x1, x2, y1, y2);
    giza_save();
    giza_set_colour_index(GPLK_GRID);
    giza_set_line_width(0.6);
    giza_set_line_style(GIZA_LS_DOT);
    giza_box_float("BG", 0, 0, "BG", 0, 0);
    if (y1 < 0 && y2 > 0) {
        giza_set_colour_index(GPLK_MUTED);
        giza_set_line_width(0.7);
        giza_move_float(x1, 0);
        giza_draw_float(x2, 0);
    }
    giza_restore();
}

// Callers use a 0..1 UI window. Hit boxes live separately in checkMenu3;
// moving a button here alone does not move its clickable area.
static void gplkDrawButton(float x1, float x2, float y1, float y2,
                           const char *label, bool active = false, bool enabled = true,bool literal=false,float fontScale=1)
{
    giza_save();
    giza_set_fill(GIZA_FILL_SOLID);
    giza_set_colour_index(active ? GPLK_ACCENT : GPLK_GRID);
    giza_rectangle_rounded_float(x1, x2, y1, y2, (y2-y1)*0.15f);
    giza_set_colour_index(!enabled ? GPLK_MUTED : active ? 0 : 1);
    float vx1,vx2,vy1,vy2,wx1,wx2,wy1,wy2,w,h;
    giza_get_viewport_float(GIZA_UNITS_NORMALIZED,&vx1,&vx2,&vy1,&vy2);
    giza_get_window_float(&wx1,&wx2,&wy1,&wy2);
    giza_get_paper_size_float(GIZA_UNITS_PIXELS,&w,&h);
    const float pixelsX=(x2-x1)/(wx2-wx1)*(vx2-vx1)*w;
    const float pixelsY=(y2-y1)/(wy2-wy1)*(vy2-vy1)*h;
    const float fit=std::min(1.0f,std::min(pixelsX/std::max(1.0f,float(std::strlen(label))*6.8f),pixelsY/15));
    giza_set_character_height_float(gplkUiCharacterHeight(0.85f)*fit*fontScale);
    if (literal) {
        const float width=gplkLiteralText(0,0,label,false);
        gplkLiteralText((x1+x2-width)*.5f,y1+(y2-y1)*.28f,label);
    } else giza_ptext_float((x1+x2)*0.5f, y1+(y2-y1)*0.28f, 0, 0.5f, label);
    giza_restore();
}

/**
 * @brief Draw residual uncertainty bars in the current series colour.
 * @param n Number of samples in each borrowed array.
 * @param x World-coordinate positions.
 * @param lo Absolute lower y bounds, not error magnitudes.
 * @param hi Absolute upper y bounds, in the same units as the current y window.
 * @param caps Native Giza terminal-cap size; zero omits caps.
 */
static void gplkDrawErrorBars(int n, const float *x, const float *lo,
                             const float *hi, float caps)
{
    if (gplkTheme == GPLK_CLASSIC) {
        giza_error_bars_vert_float(n, x, lo, hi, caps);
        return;
    }
    giza_save();
    int colour;
    double r, g, b;
    giza_get_colour_index(&colour);
    giza_get_colour_representation(colour, &r, &g, &b);
    giza_set_colour_representation_alpha(GPLK_ERROR, r, g, b, 0.42);
    giza_set_colour_index(GPLK_ERROR);
    giza_set_line_width(0.65);
    giza_error_bars_vert_float(n, x, lo, hi, caps);
    giza_restore();
}

static void gplkDrawFrequencyLegend(int n, const float *lo, const float *hi, const int *colours)
{
    if (gplkTheme==GPLK_CLASSIC || n==0) return;
    giza_save();
    giza_set_character_height(0.6);
    giza_set_colour_index(GPLK_MUTED);
    giza_text(0.62,0.74,"Frequency (MHz)");
    for (int i=0; i<n && i<5; ++i) {
        char label[48];
        if (hi[i]>=100000) std::snprintf(label,sizeof(label),"%.3gk+",lo[i]/1000);
        else if (hi[i]>=1000) std::snprintf(label,sizeof(label),"%.3g-%.3gk",lo[i]/1000,hi[i]/1000);
        else if (lo[i]<=0) std::snprintf(label,sizeof(label),"<%.0f",hi[i]);
        else std::snprintf(label,sizeof(label),"%.0f-%.0f",lo[i],hi[i]);
        giza_set_colour_index(colours[i]);
        giza_text(0.62+0.076*i,0.58,label);
    }
    giza_restore();
}

/**
 * @brief Open/select a device and initialize its colours, theme and panel grid.
 * @return True on success; false after reporting an open failure.
 * Does not close an existing device. The caller owns closing/reselecting it;
 * a child spectrum must not close all devices when returning to the main plot.
 */
static bool gplkOpenDevice(const char *device, int columns, int rows)
{
    if (giza_open_device(device, "gplk") <= 0) {
        std::fprintf(stderr, "gplk: cannot open Giza device '%s'\n", device);
        return false;
    }
    giza_set_colour_palette(GIZA_COLOUR_PALETTE_PGPLOT);
    char type[32] = {0};
    int length = sizeof(type);
    giza_query_device("type", type, &length);
    const bool screen = giza_device_has_cursor() || std::strcmp(type, "/png") == 0 ||
                        std::strcmp(type, "/null") == 0;
    const double background = screen ? 0.0 : 1.0;
    giza_set_colour_representation(0, background, background, background);
    giza_set_colour_representation(1, 1.0-background, 1.0-background, 1.0-background);
    gplkApplyTheme();
    giza_stop_prompting();
    if (columns != 1 || rows != 1) {
        giza_subpanel(columns, rows);
        giza_set_panel(1, 1);
    }
    return true;
}

/** Keep device width; set height/width when aspect>0, otherwise retain its size. */
static void gplkSetAspect(float aspect)
{
    float width, height;
    giza_get_paper_size_float(GIZA_UNITS_INCHES, &width, &height);
    if (!(width > 0.0f)) width = 8.0f;
    if (aspect > 0.0f) height = width * aspect;
    if (!(height > 0.0f)) height = width * 0.8f;
    giza_set_paper_size_float(GIZA_UNITS_INCHES, width, height);
}

static void gplkSetFont(int font)
{
    switch (font) {
    case 2: giza_set_font("serif"); break;
    case 3: giza_set_font_italic("serif"); break;
    case 4: giza_set_font("cursive"); break;
    default: giza_set_font("sans-serif"); break;
    }
}

static std::string gplkColourName(const char *name)
{
    std::string result;
    for (; *name; ++name) {
        unsigned char c = static_cast<unsigned char>(*name);
        if (!std::isspace(c)) result += static_cast<char>(std::tolower(c));
    }
    return result;
}

// Giza sets RGB values directly. Resolve setup-file names without opening an
// X display, so publication plots also work on machines without a display.
static bool gplkSetNamedColour(int index, const char *name)
{
    const std::string key = gplkColourName(name);
    struct Colour { const char *name; unsigned int rgb; };
    static const Colour colours[] = {
        {"black", 0x000000}, {"white", 0xffffff}, {"red", 0xff0000},
        {"green", 0x00ff00}, {"blue", 0x0000ff}, {"cyan", 0x00ffff},
        {"magenta", 0xff00ff}, {"yellow", 0xffff00}, {"orange", 0xffa500},
        {"gray", 0x808080}, {"grey", 0x808080}, {"darkgray", 0xa9a9a9},
        {"darkgrey", 0xa9a9a9}, {"lightgray", 0xd3d3d3},
        {"lightgrey", 0xd3d3d3}, {"brown", 0xa52a2a}, {"pink", 0xffc0cb},
        {"purple", 0xa020f0}, {"darkgreen", 0x006400}, {"skyblue", 0x87ceeb}
    };
    unsigned int rgb = 0;
    bool found = false;
    if (key.size() == 7 && key[0] == '#' &&
        key.find_first_not_of("0123456789abcdef", 1) == std::string::npos) {
        rgb = static_cast<unsigned int>(std::strtoul(key.c_str() + 1, NULL, 16));
        found = true;
    }
    for (unsigned int i = 0; !found && i < sizeof(colours)/sizeof(colours[0]); ++i) {
        if (key == colours[i].name) { rgb = colours[i].rgb; found = true; }
    }
    if (!found) {
        const char *paths[] = {"/usr/share/X11/rgb.txt", "/usr/lib/X11/rgb.txt"};
        for (unsigned int i = 0; !found && i < sizeof(paths)/sizeof(paths[0]); ++i) {
            FILE *file = std::fopen(paths[i], "r");
            if (!file) continue;
            char line[512], colour[512];
            unsigned int r, g, b;
            while (std::fgets(line, sizeof(line), file)) {
                if (std::sscanf(line, "%u %u %u %511[^\n]", &r, &g, &b, colour) == 4 &&
                    r <= 255 && g <= 255 && b <= 255 && key == gplkColourName(colour)) {
                    rgb = (r << 16) | (g << 8) | b;
                    found = true;
                    break;
                }
            }
            std::fclose(file);
        }
    }
    if (!found) {
        std::fprintf(stderr, "gplk: unknown colour '%s'; keeping colour index %d\n", name, index);
        return false;
    }
    giza_set_colour_representation(index, ((rgb >> 16) & 255) / 255.0,
                                   ((rgb >> 8) & 255) / 255.0, (rgb & 255) / 255.0);
    return true;
}

/** Return the full surface extent in native GIZA_UNITS_* units, origin at zero. */
static void gplkGetSurfaceSize(int units, float *x1, float *x2, float *y1, float *y2)
{
    *x1 = *y1 = 0.0f;
    giza_get_paper_size_float(units, x2, y2);
}

static void gplkGetViewportPixels(float *x1, float *x2, float *y1, float *y2)
{
    // Some Giza versions reflect physical viewport queries about the top of
    // the surface. gplk hit testing needs a bottom-left origin; derive it from
    // the normalised viewport instead of inheriting that reflection.
    float width, height;
    giza_get_paper_size_float(GIZA_UNITS_PIXELS, &width, &height);
    giza_get_viewport_float(GIZA_UNITS_NORMALIZED, x1, x2, y1, y2);
    *x1 *= width;
    *x2 *= width;
    *y1 *= height;
    *y2 *= height;
}

// Overlay colours are semantic roles, independent of the residuals' palette.
static void gplkPrepareOverlayColours()
{
    double r,g,b;
    giza_get_colour_representation(0,&r,&g,&b);
    const bool dark=0.2126*r+0.7152*g+0.0722*b<0.5;
    gplkSetRGB(GPLK_BINS,dark ? 0xffdf60 : 0x9b4400);
    gplkSetRGB(GPLK_MODEL,dark ? 0xfa79e1 : 0x962580);
    gplkSetRGB(GPLK_SELECTED,dark ? 0xffffff : 0x101820);
    gplkSetRGB(GPLK_SOLAR,dark ? 0x9ce76c : 0x456b00);
    gplkSetRGB(GPLK_NOISE_OVERLAY,dark ? 0x75ddff : 0x005f87);
}

/**
 * @brief Draw binned means above residuals/model curves, with opaque error bars.
 * @param n Number of samples; a nonpositive count is a no-op.
 * @param x Bin centres in current x-world units.
 * @param y Bin means in current y-world units.
 * @param error Symmetric y-error magnitudes, in y-world units (unlike lo/hi above).
 * Marker radius starts in pixels and is converted to world offsets, so zooming
 * does not shrink the overlay. This helper styles existing bin statistics only.
 */
static void gplkDrawBinnedOverlay(int n, const float *x, const float *y, const float *error)
{
    if (n<=0) return;
    gplkPrepareOverlayColours();
    giza_save();
    giza_set_line_style(GIZA_LS_SOLID);
    float x1,x2,y1,y2,vx1,vx2,vy1,vy2;
    giza_get_window_float(&x1,&x2,&y1,&y2);
    gplkGetViewportPixels(&vx1,&vx2,&vy1,&vy2);
    float height;
    giza_get_character_height_float(&height);
    const float radius=6.0f*std::max(1.0f,std::min(2.0f,height));
    const float rx=(x2-x1)*radius/std::max(1.0f,std::fabs(vx2-vx1));
    const float ry=(y2-y1)*radius/std::max(1.0f,std::fabs(vy2-vy1));
    const int halo=gplkTheme==GPLK_CLASSIC ? 0 : GPLK_PANEL;
    // Knock out a small halo before drawing opaque error bars and diamonds.
    giza_set_colour_index(halo);
    giza_set_line_width(5.2);
    giza_error_bars_float(6,n,x,y,error,1.3f);
    giza_set_colour_index(GPLK_BINS);
    giza_set_line_width(2.2);
    giza_error_bars_float(6,n,x,y,error,1.3f);
    for (int i=0; i<n; ++i) {
        const float outerX[]={x[i],x[i]+1.3f*rx,x[i],x[i]-1.3f*rx};
        const float outerY[]={y[i]+1.3f*ry,y[i],y[i]-1.3f*ry,y[i]};
        const float innerX[]={x[i],x[i]+rx,x[i],x[i]-rx};
        const float innerY[]={y[i]+ry,y[i],y[i]-ry,y[i]};
        giza_set_fill(GIZA_FILL_SOLID);
        giza_set_colour_index(halo);
        giza_polygon_float(4,outerX,outerY);
        giza_set_colour_index(GPLK_BINS);
        giza_polygon_float(4,innerX,innerY);
        giza_set_fill(GIZA_FILL_HOLLOW);
        giza_set_colour_index(GPLK_SELECTED);
        giza_set_line_width(1.1);
        giza_polygon_float(4,innerX,innerY);
    }
    giza_restore();
}

/** Draw a haloed model polyline in input order; arrays are in current world units. */
static void gplkDrawOverlayLine(int n, const float *x, const float *y, int style, int colour=GPLK_MODEL)
{
    if (n<=0) return;
    gplkPrepareOverlayColours();
    giza_save();
    giza_set_line_style(style);
    giza_set_colour_index(gplkTheme==GPLK_CLASSIC ? 0 : GPLK_PANEL);
    giza_set_line_width(5.0);
    giza_line_float(n,x,y);
    giza_set_colour_index(colour);
    giza_set_line_width(2.3);
    giza_line_float(n,x,y);
    giza_restore();
}

/** Preserve the selection group's Giza symbol; emphasize its outline and halo. */
static void gplkDrawHighlightedOverlay(int n, const float *x, const float *y, int symbol)
{
    if (n<=0) return;
    gplkPrepareOverlayColours();
    giza_save();
    giza_set_character_height(gplkPlotCharacterHeight(1.15f));
    giza_set_colour_index(gplkTheme==GPLK_CLASSIC ? 0 : GPLK_PANEL);
    giza_set_line_width(5.0);
    gplkDrawPoints(n,x,y,symbol);
    giza_set_colour_index(GPLK_SELECTED);
    giza_set_line_width(2.0);
    gplkDrawPoints(n,x,y,symbol);
    giza_restore();
}

static void gplkDrawOverlayLegend(bool bins, bool shapiro, bool fitwaves, bool selected,
                                  bool solar=false, const char *noise=NULL)
{
    std::string text;
    if (bins) text="Bins: diamonds";
    if (shapiro) text+=(text.empty() ? "" : " | ")+std::string("Shapiro: dashed");
    if (fitwaves) text+=(text.empty() ? "" : " | ")+std::string("FITWAVES: dashed");
    if (selected) text+=(text.empty() ? "" : " | ")+std::string("Selected: outlined");
    if (solar) text+=(text.empty() ? "" : " | ")+std::string("Solar: green");
    if (noise) text+=(text.empty() ? "" : " | ")+std::string(noise)+": squares";
    if (text.empty()) return;
    giza_save();
    giza_set_character_height(gplkUiCharacterHeight(text.size()>80 ? 0.52f : 0.62f));
    giza_set_colour_index(1);
    giza_annotate("T",0.45,0,0,text.c_str());
    giza_restore();
}

// Giza returns zero on success. Use the double entry point so failed input
// never copies uninitialised float coordinates from the library's wrappers.
static bool gplkReadBand(int mode, int move, float anchorX, float anchorY,
                        float *x, float *y, char *key)
{
    double dx = anchorX, dy = anchorY;
    *key = 'q';
    if (!giza_device_has_cursor()) return false;
    // end_buffer already presents the frame; flushing here would present it
    // a second time for every input event.
    if (giza_band(mode, move, anchorX, anchorY, &dx, &dy, key) != 0) {
        *key = 'q';
        return false;
    }
    *x = static_cast<float>(dx);
    *y = static_cast<float>(dy);
    return true;
}

static void gplkResizeMotion(double *,double *,int *) {}

/** Main idle input only: newer Giza X11 drivers opt into resize notifications
 * through a motion callback and return synthetic 'r'. Never dispatch that as
 * gplk's destructive restart/refit command. Nested rubber-band input keeps its
 * ordinary semantics. Old drivers can use R to redraw at their new size. */
static bool gplkReadEvent(float *x,float *y,char *key)
{
    float beforeW,beforeH,afterW,afterH;
    giza_get_paper_size_float(GIZA_UNITS_PIXELS,&beforeW,&beforeH);
    giza_set_motion_callback(gplkResizeMotion);
    const bool ok=gplkReadBand(GIZA_BAND_NONE,0,0,0,x,y,key);
    giza_end_motion_callback();
    giza_get_paper_size_float(GIZA_UNITS_PIXELS,&afterW,&afterH);
    if (ok && *key!='q' && (std::fabs(afterW-beforeW)>0.5f || std::fabs(afterH-beforeH)>0.5f))
        *key=0; // redraw only; ignore coordinates from the obsolete viewport
    return ok;
}

#endif
