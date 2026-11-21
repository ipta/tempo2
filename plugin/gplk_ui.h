/** Native Giza responsive controls. No timing-model or fit ownership here.
 * Menu state is local to doPlot; only an explicit jump click changes fit flags.
 * The colour-group cache must be rebuilt after replacing observations.
 * Distributed under TEMPO2's GPLv3-or-later licence. */
#ifndef TEMPO2_GPLK_UI_H
#define TEMPO2_GPLK_UI_H
#include <map>
#include <vector>
#include <sstream>
#include "gplk_jumps.h"

struct GplkColourGroups {
    std::string mode;
    std::vector<std::string> choices,values;
    std::vector<int> group;
    std::vector<std::vector<int> > plotGroups;
    bool valid;
    GplkColourGroups():mode("@frequency"),valid(false) {}
    bool categorical() const { return mode!="@frequency" && mode!="@setup" && mode!="@single"; }
    static int colour(int index) { return 2+index%11; }
    void groupPositions(int count,const int *ids) {
        plotGroups.resize(categorical() ? values.size() : 1);
        for (unsigned int i=0;i<plotGroups.size();++i) plotGroups[i].clear();
        for (int i=0;i<count;++i) plotGroups[categorical() ? group[ids[i]] : 0].push_back(i);
    }
    void rebuild(const pulsar &p) {
        choices.clear(); choices.push_back("@frequency"); choices.push_back("@single");
        choices.push_back("@setup"); choices.push_back("@site");
        std::map<std::string,int> flags,groups;
        values.clear(); group.assign(p.nobs,0);
        for (int i=0;i<p.nobs;++i) {
            std::string value=mode=="@site" ? p.obsn[i].telID : "(missing)";
            for (int k=0;k<p.obsn[i].nFlags;++k) {
                const std::string flag=p.obsn[i].flagID[k];
                if (flag!="-selectAll") flags[flag]=1;
                if (flag==mode) value=p.obsn[i].flagVal[k];
            }
            if (categorical()) {
                std::map<std::string,int>::iterator found=groups.find(value);
                if (found==groups.end()) {
                    group[i]=values.size(); groups[value]=group[i]; values.push_back(value);
                } else group[i]=found->second;
            }
        }
        for (std::map<std::string,int>::const_iterator f=flags.begin();f!=flags.end();++f)
            choices.push_back(f->first);
        valid=true;
    }
    static std::string label(const std::string &m) {
        if (m=="@frequency") return "Frequency";
        if (m=="@single") return "Single color";
        if (m=="@setup") return "Setup rules";
        if (m=="@site") return "Observatory";
        return m;
    }
};

struct GplkFooterItem {
    std::string label;
    int kind,index,colour;
    bool active,fixed;
    GplkFooterItem(const std::string &s,int k,int i,int c,bool a,bool f=false):label(s),kind(k),index(i),colour(c),active(a),fixed(f) {}
};
struct GplkUiState {
    std::vector<GplkFooterItem> jumps,legend;
    int page,legendPage,dropdownPage;
    bool dropdown,initialLayout;
    GplkUiState():page(0),legendPage(0),dropdownPage(0),dropdown(false),initialLayout(true) {}
    void rebuild(const pulsar &p,const GplkColourGroups &colours) {
        jumps.clear(); legend.clear();
        for (int kind=1;kind<=2;++kind) {
            const int n=kind==1 ? p.nJumps : p.nfdJumps;
            for (int i=1;i<=n;++i) {
                const char *raw=kind==1 ? p.jumpStr[i] : p.fdjumpStr[i];
                const GplkJumpDescription description=gplkDescribeJump(raw);
                if (description.match=="-1") continue;
                const bool fixed=kind==1 && description.correction;
                const bool active=!fixed && (kind==1 ? p.fitJump[i]!=0 : p.fitfdJump[i]!=0);
                std::string label=description.label;
                if (kind==2) {
                    std::ostringstream prefix;
                    if (p.fdjumpIdx[i]==-2) prefix<<"FD DM "; else prefix<<"FD"<<p.fdjumpIdx[i]<<" ";
                    label=prefix.str()+label;
                }
                jumps.push_back(GplkFooterItem(label,kind,i,active ? GPLK_ACCENT : GPLK_MUTED,active,fixed));
            }
        }
        for (unsigned int i=0;i<colours.values.size();++i)
            legend.push_back(GplkFooterItem(colours.values[i],0,i,GplkColourGroups::colour(i),false));
        if (colours.mode=="@single")
            legend.push_back(GplkFooterItem("All observations",0,0,GplkColourGroups::colour(0),false));
    }
};

// Fit long flag values without changing their stored value or click identity.
static std::string gplkElide(const std::string &s,int characters) {
    return int(s.size())<=characters ? s : s.substr(0,std::max(1,characters-3))+"...";
}

static void gplkUpdateLayout(int entries,bool grow,int legendEntries=0,bool initial=false)
{
    float w,h;
    giza_get_paper_size_float(GIZA_UNITS_PIXELS,&w,&h);
    gplkLayout.update(w,h,entries);
    if (grow && giza_device_has_cursor()) {
        const float requestedWidth=initial ? std::max(w,1050.0f) : w;
        gplkLayout.update(requestedWidth,h,entries);
        const int rows=(entries+gplkLayout.columns-1)/gplkLayout.columns;
        // Grow once when the list changes, never fight a user's resize drag.
        const float wanted=std::min(1000.0f,std::max(620.0f,std::max((82+24*rows)/0.42f,
                                      112+82+24*rows+std::max(376.0f,116.0f+24*legendEntries))));
        if (wanted>h+1 || requestedWidth>w+1) {
            giza_set_paper_size_float(GIZA_UNITS_PIXELS,requestedWidth,std::max(h,wanted));
            giza_change_page();
            giza_get_paper_size_float(GIZA_UNITS_PIXELS,&w,&h);
            gplkLayout.update(w,h,entries);
        }
    }
}

static int gplkFooterCapacity() { return gplkLayout.columns*gplkLayout.visibleRows; }
static GplkRect gplkFooterCell(int slot) {
    const float row=24/gplkLayout.height/gplkLayout.bottom.top;
    const float top=1-82/gplkLayout.height/gplkLayout.bottom.top;
    const int col=slot%gplkLayout.columns,r=slot/gplkLayout.columns;
    return GplkRect(float(col)/gplkLayout.columns+0.004f,float(col+1)/gplkLayout.columns-0.004f,
                   top-(r+1)*row,top-r*row-0.012f);
}

static void gplkDrawFooter(const pulsar &p,int yplot,GplkUiState &ui)
{
    gplkLayout.bottom.select(); giza_set_window_float(0,1,0,1); gplkDrawPanel(0,1,0,1);
    const float row=24/gplkLayout.height/gplkLayout.bottom.top;
    const char *labels[]={"Quit","Clear","Measure","Help","Red","DM","Chrom","Total","Clean","Residual"};
    const int modes[]={GPLK_NOISE_RED,GPLK_NOISE_DM,GPLK_NOISE_CHROM,GPLK_NOISE_TOTAL,GPLK_NOISE_CLEAN,2};
    for (int b=0;b<10;++b)
        gplkDrawButton(b*.1f+.005f,b*.1f+.085f,1-row+.006f,.995f,labels[b],b>=4 && yplot==modes[b-4],
                      b<4 || b==9 || gplkNoiseAvailable(p,modes[b-4]));
    const char *palettes[]={"Ocean","Contrast","Sunset"};
    const char *second[]={"List","Print","Highlight",gplkTheme==GPLK_DARK ? "Light" : "Dark",palettes[gplkPalette],"Spectrum"};
    for (int b=0;b<6;++b)
        gplkDrawButton(b*.1f+.005f,b*.1f+.085f,1-2*row+.006f,1-row-.008f,second[b]);
    const int capacity=gplkFooterCapacity();
    const int pages=capacity>0 ? std::max(1,(int(ui.jumps.size())+capacity-1)/capacity) : 1;
    ui.page=std::min(ui.page,pages-1);
    // A separate labelled subarea below the toolbar owns jumps. It has
    // its own background, divider, wrapping and paging; toolbar height is fixed.
    const float jumpsTop=gplkLayout.jumps.top/gplkLayout.bottom.top;
    giza_save();
    giza_set_colour_index(0); giza_set_fill(GIZA_FILL_SOLID);
    giza_rectangle_float(0,1,0,jumpsTop);
    giza_set_colour_index(GPLK_GRID);
    giza_move_float(0,jumpsTop); giza_draw_float(1,jumpsTop);
    giza_set_colour_index(GPLK_MUTED);
    giza_set_character_height_float(gplkUiCharacterHeight(.85f));
    char caption[160];
    std::snprintf(caption,sizeof(caption),"Jumps | %u entries | page %d/%d",unsigned(ui.jumps.size()),ui.page+1,pages);
    int fixedCount=0;
    for (unsigned int i=0;i<ui.jumps.size();++i) if (ui.jumps[i].fixed) ++fixedCount;
    if (fixedCount)
        std::snprintf(caption,sizeof(caption),"Jumps | %u entries | page %d/%d | %d fixed (applied)",
                      unsigned(ui.jumps.size()),ui.page+1,pages,fixedCount);
    if (capacity==0 && !ui.jumps.empty())
        std::snprintf(caption,sizeof(caption),"Jumps | %u entries | enlarge window to view",unsigned(ui.jumps.size()));
    giza_text_float(.008f,jumpsTop-18/gplkLayout.height/gplkLayout.bottom.top,caption);
    giza_restore();
    if (pages>1) {
        gplkDrawButton(.865f,.92f,jumpsTop-row+.006f,jumpsTop-.008f,"<",false,ui.page>0);
        gplkDrawButton(.925f,.985f,jumpsTop-row+.006f,jumpsTop-.008f,">",false,ui.page+1<pages);
    }
    for (int slot=0;slot<capacity && ui.page*capacity+slot<int(ui.jumps.size());++slot) {
        const GplkFooterItem &item=ui.jumps[ui.page*capacity+slot];
        const GplkRect cell=gplkFooterCell(slot);
        giza_save();
        if (item.fixed) {
            // Read-only corrections are text only: no fill, border, badge box,
            // or button chrome. Keep the full flag name readable beside "fixed".
            giza_set_colour_index(GPLK_MUTED);
            giza_set_character_height_float(gplkUiCharacterHeight(.70f));
            gplkLiteralText(cell.left+.006f,cell.bottom+(cell.top-cell.bottom)*.28f,"fixed");
        } else gplkDrawButton(cell.left,cell.right,cell.bottom,cell.top,"",item.active);
        giza_set_colour_index(item.active ? 0 : 1);
        giza_set_character_height_float(gplkUiCharacterHeight(.85f));
        const float inset=item.fixed ? 42/gplkLayout.width : .006f;
        const int nch=std::max(6,int((cell.right-cell.left-inset)*gplkLayout.width/7)-1);
        const std::string label=gplkElide(item.label,nch);
        gplkLiteralText(cell.left+inset,cell.bottom+(cell.top-cell.bottom)*.28f,label);
        giza_restore();
    }
}

// Right-panel geometry is in full-surface coordinates for both rendering and
// input. The dropdown opens downward over its legend, never over jump buttons.
static GplkRect gplkColourSelector() {
    const GplkRect &panel=gplkLayout.colours;
    return GplkRect(panel.left+8/gplkLayout.width,panel.right-8/gplkLayout.width,
                   panel.top-56/gplkLayout.height,panel.top-30/gplkLayout.height);
}
static int gplkLegendCapacity() {
    return std::max(0,int(((gplkLayout.colours.top-gplkLayout.colours.bottom)*gplkLayout.height-116)/24));
}
static GplkRect gplkLegendPager(bool next) {
    const GplkRect &panel=gplkLayout.colours;
    const float mid=(panel.left+panel.right)*.5f;
    return GplkRect(next ? mid+4/gplkLayout.width : panel.left+8/gplkLayout.width,
                   next ? panel.right-8/gplkLayout.width : mid-4/gplkLayout.width,
                   panel.bottom+6/gplkLayout.height,panel.bottom+30/gplkLayout.height);
}
static void gplkDrawColourPanel(const GplkColourGroups &colours,GplkUiState &ui) {
    giza_save(); giza_set_viewport_float(0,1,0,1); giza_set_window_float(0,1,0,1);
    const GplkRect &panel=gplkLayout.colours;
    gplkDrawPanel(panel.left,panel.right,panel.bottom,panel.top);
    giza_set_colour_index(1); giza_set_character_height_float(gplkUiCharacterHeight(.95f));
    giza_text_float(panel.left+10/gplkLayout.width,panel.top-20/gplkLayout.height,"Color by");
    const GplkRect selector=gplkColourSelector();
    const std::string label=gplkElide(GplkColourGroups::label(colours.mode),30)+" v";
    gplkDrawButton(selector.left,selector.right,selector.bottom,selector.top,label.c_str(),ui.dropdown,true,true);
    const int capacity=gplkLegendCapacity();
    const int pages=capacity>0 ? std::max(1,(int(ui.legend.size())+capacity-1)/capacity) : 1;
    ui.legendPage=std::min(ui.legendPage,pages-1);
    char caption[96];
    std::snprintf(caption,sizeof(caption),"%u groups | %d/%d",unsigned(ui.legend.size()),ui.legendPage+1,pages);
    giza_set_colour_index(GPLK_MUTED); giza_set_character_height_float(gplkUiCharacterHeight(.80f));
    giza_text_float(panel.left+10/gplkLayout.width,panel.top-78/gplkLayout.height,caption);
    for (int slot=0;slot<capacity && ui.legendPage*capacity+slot<int(ui.legend.size());++slot) {
        const GplkFooterItem &item=ui.legend[ui.legendPage*capacity+slot];
        const float top=panel.top-(88+24*slot)/gplkLayout.height;
        giza_set_colour_index(item.colour); giza_set_fill(GIZA_FILL_SOLID);
        giza_rectangle_float(panel.left+10/gplkLayout.width,panel.left+20/gplkLayout.width,
                             top-16/gplkLayout.height,top-6/gplkLayout.height);
        giza_set_colour_index(1); giza_set_character_height_float(gplkUiCharacterHeight(.85f));
        const int nch=std::max(5,int((panel.right-panel.left)*gplkLayout.width-44)/7);
        gplkLiteralText(panel.left+30/gplkLayout.width,top-16/gplkLayout.height,gplkElide(item.label,nch));
    }
    if (pages>1) {
        const GplkRect prev=gplkLegendPager(false),next=gplkLegendPager(true);
        gplkDrawButton(prev.left,prev.right,prev.bottom,prev.top,"<",false,ui.legendPage>0);
        gplkDrawButton(next.left,next.right,next.bottom,next.top,">",false,ui.legendPage+1<pages);
    }
    giza_restore();
}
static GplkRect gplkDropdownRect(int count) {
    const GplkRect selector=gplkColourSelector();
    const float top=selector.bottom-4/gplkLayout.height;
    return GplkRect(selector.left,selector.right,top-24*count/gplkLayout.height,top);
}
static int gplkDropdownRows() {
    return std::max(2,std::min(16,int((gplkColourSelector().bottom-gplkLayout.colours.bottom)*gplkLayout.height/24)-1));
}
static void gplkDrawColourDropdown(const GplkColourGroups &colours,GplkUiState &ui) {
    if (!ui.dropdown) return;
    giza_save(); giza_set_viewport_float(0,1,0,1); giza_set_window_float(0,1,0,1);
    const int rows=gplkDropdownRows();
    ui.dropdownPage=std::min(ui.dropdownPage,(int(colours.choices.size())-1)/rows);
    const int n=std::min(rows,int(colours.choices.size())-ui.dropdownPage*rows);
    const bool paged=int(colours.choices.size())>rows;
    const GplkRect box=gplkDropdownRect(n+(paged ? 1 : 0));
    const float row=24/gplkLayout.height;
    gplkDrawPanel(box.left,box.right,box.bottom,box.top);
    for (int i=0;i<n;++i) {
        const std::string &mode=colours.choices[ui.dropdownPage*rows+i];
        const std::string label=gplkElide(GplkColourGroups::label(mode),32);
        gplkDrawButton(box.left+.003f,box.right-.003f,box.top-(i+1)*row+.002f,box.top-i*row-.002f,
                      label.c_str(),mode==colours.mode,true,true);
    }
    if (paged) {
        gplkDrawButton(box.left+.003f,(box.left+box.right)/2-.003f,box.bottom+.002f,box.bottom+row-.002f,"< Previous");
        gplkDrawButton((box.left+box.right)/2+.003f,box.right-.003f,box.bottom+.002f,box.bottom+row-.002f,"Next >");
    }
    giza_restore(); // caller restores plot viewport/window explicitly
}

// Returns true when the modern footer/dropdown consumed the event. Other
// panels continue through legacy command dispatch using the shared layout.
static bool gplkHandleUi(float nx,float ny,char &key,pulsar &p,GplkColourGroups &colours,GplkUiState &ui,bool &grow)
{
    if (ui.dropdown) {
        if (key==27 || key=='q') { ui.dropdown=false; key=0; return true; }
        if (key!='A') { key=0; return true; }
        const int rows=gplkDropdownRows();
        ui.dropdownPage=std::min(ui.dropdownPage,(int(colours.choices.size())-1)/rows);
        const int n=std::min(rows,int(colours.choices.size())-ui.dropdownPage*rows);
        const bool paged=int(colours.choices.size())>rows;
        const GplkRect box=gplkDropdownRect(n+(paged ? 1 : 0));
        if (box.contains(nx,ny)) {
            const int row=int((box.top-ny)*gplkLayout.height/24);
            if (row<n) {
                colours.mode=colours.choices[ui.dropdownPage*rows+row]; colours.rebuild(p);
                ui.dropdown=false; ui.legendPage=0; grow=true;
            } else if (paged) {
                const int pages=(colours.choices.size()+rows-1)/rows;
                ui.dropdownPage=std::max(0,std::min(pages-1,ui.dropdownPage+(nx<(box.left+box.right)/2 ? -1 : 1)));
            }
        } else ui.dropdown=false;
        key=0; return true;
    }
    if (key=='A' && gplkLayout.colours.contains(nx,ny)) {
        if (gplkColourSelector().contains(nx,ny)) { ui.dropdown=true; ui.dropdownPage=0; }
        else if (gplkLegendPager(false).contains(nx,ny)) ui.legendPage=std::max(0,ui.legendPage-1);
        else if (gplkLegendPager(true).contains(nx,ny) && gplkLegendCapacity()>0 &&
                 (ui.legendPage+1)*gplkLegendCapacity()<int(ui.legend.size())) ++ui.legendPage;
        key=0; return true;
    }
    if (key!='A' || !gplkLayout.bottom.contains(nx,ny)) return false;
    const float lx=gplkLayout.bottom.localX(nx),ly=gplkLayout.bottom.localY(ny);
    const int row=int((1-ly)*gplkLayout.bottom.top*gplkLayout.height/24),col=int(lx*10);
    if (row==0) {
        if (col==0) key='q';
        if (col==3) key='h';
        // Other toolbar actions retain checkMenu3 dispatch in the shared layout.
        return false;
    }
    if (row==1 && lx>=.6f) {
        key=0; return true;
    }
    if (row<2) return false;
    const float jumpsTop=gplkLayout.jumps.top/gplkLayout.bottom.top;
    if (lx>=.865f && ly<jumpsTop && ly>=jumpsTop-24/gplkLayout.height/gplkLayout.bottom.top) {
        if (lx<.925f) ui.page=std::max(0,ui.page-1);
        else if (gplkFooterCapacity()>0 && (ui.page+1)*gplkFooterCapacity()<int(ui.jumps.size())) ++ui.page;
        key=0; return true;
    }
    const int capacity=gplkFooterCapacity();
    for (int slot=0;slot<capacity && ui.page*capacity+slot<int(ui.jumps.size());++slot) {
        if (!gplkFooterCell(slot).contains(lx,ly)) continue;
        const GplkFooterItem &item=ui.jumps[ui.page*capacity+slot];
        if (item.fixed) break; // text entries have no button action
        if (item.kind==1) gplkToggleJumpFit(p,item.index);
        if (item.kind==2) { p.fitfdJump[item.index]=!p.fitfdJump[item.index]; p.fdjumpValErr[item.index]=0; }
        break;
    }
    key=0; return true;
}
#endif
