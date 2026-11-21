Git INSTALLATION README
=======================

0. Contents
	1. What this package is    
	2. Quick Guide    
	3. Requirements    
	4. Detailed instalation guide    
	5. Plugins    
	6. Changes from old makefile    
	7. Installation troubleshooting    
	8. Bugs and feature requests


1. What this package is
------------------------
You (or someone else) have checked out tempo2 from the Git
(https://bitbucket.org/psrsoft/tempo2)

This is the best way to get the latest/cutting edge version, and develop your
own additions to the tempo2 code or via plugins.

For more information on tempo2 see:
http://www.atnf.csiro.au/research/pulsar/tempo2/

This requires the gnu autotools. If you don't have or don't want to install 
autotools, we recommend you install the latest distributed release from
http://www.atnf.csiro.au/research/pulsar/tempo2/
or use PSRSOFT to install tempo2:
http://www.pulsarastronomy.net/wiki/Software/PSRSoft

2. Quick Guide
---------------
Bootstrap the build system:

    ./bootstrap

setup the tempo2 runtime dir

    cp -r T2runtime /usr/share/tempo2/
    export TEMPO2=/usr/share/tempo2/

Configure:

    ./configure [[--prefix=/your/install/path]]

use --prefix to set the path you want to install the binaries and libraries

Make and install...

    make && make install

You will probably want to build the default plugins (gplk, etc). Do this with:

    make plugins && make plugins-install


And you're done.

3. Requirements
---------------
Tempo2 requires the following:

 - A fortran 77 compiler (tested with gfortran).
 - A C compiler (tested with gcc).

The `gplk` plugin requires native [Giza](https://danieljprice.github.io/giza/)
(`giza.h` and `libgiza`). It does not use Giza's PGPLOT compatibility libraries.
Other plotting plugins still require PGPLOT.

5. Plugins
----------
The bootstrap command will create suitible makefiles for the default set of
plugins. This is controled by the contents of the files in
./plugin/plugin_lists/

 - vanilla.plugins lists plugins to install which have no dependancies.
 - giza.plugins lists plugins to install that depend on native Giza (`gplk`).
 - pgplot.plugins lists plugins to install that are dependant on PGPLOT.
 - gsl.plugins lists plugins to install that are dependant on the GSL.

5.1 Native Giza gplk
-------------------
Build and run `gplk` in six steps:

1. Install Giza's development headers and library, with X11 support for
   interactive use. Ensure its `giza.pc` is on `PKG_CONFIG_PATH`, or note its
   installation prefix.
2. Run `./bootstrap`.
3. Set `TEMPO2` to your runtime directory and run
   `./configure --with-giza` (or `--with-giza=/path/to/giza`). Explicitly
   requesting Giza makes configuration fail if the native API cannot link.
4. Run `make`, then `make -C plugin gplk`.
5. Run `make install`, then `make -C plugin gplk-install` to install the core and
   plugin into the configured paths.
6. Run `tempo2 -gr gplk -f pulsar.par pulsar.tim`. For headless output, add
   `-grdev residuals.pdf/pdf` (also PNG, SVG and PS). Run the headless native
   renderer check with `make -C tests check TESTS=testGplkGiza`.

Giza is auto-detected by default; `--without-giza` disables `gplk`. For custom
library layouts or static linking, provide `GIZA_CFLAGS` and `GIZA_LIBS`.
Only `libgiza` is needed for `gplk`; missing PGPLOT does not disable it.
The original PGPLOT plugin is separately available as `plk`. The native Giza
implementation is named `gplk` in plugin filenames, build/install targets,
window titles, messages, helper APIs and tests. Its settings use
`gplk_setup*.dat`; `gplk-install` installs those files into the runtime
`plugin_data` directory alongside the plugin installation.

The default interactive device is `/xw`. Legacy `/xs` device strings remain
accepted by Giza. Mouse buttons, keyboard commands, setup files, refitting and
replay retain their existing gplk meanings. Setup-file `fonttype` values select
sans-serif (1), serif (2), italic serif (3), or cursive (4). Background and line
colours accept common names, `#RRGGBB`, and names from the system X11 `rgb.txt`
when available. Cursor availability is queried from the open device, so file
devices run without waiting for input.

Appearance and timing-noise views
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Developers: see the [native Giza gplk maintainer guide](documentation/GPLK_GIZA_DEVELOPER_GUIDE.md)
for source maps, extension recipes, numerical conventions and known limitations.

The native viewer defaults to a slate dark theme with a coordinated series
palette, subdued uncertainty bars, a fine grid and a zero reference line.
`-publish` and `-image` default to the light theme. Choose either explicitly:

    tempo2 -gr gplk -f pulsar.par pulsar.tim -theme light -palette ocean
    tempo2 -gr gplk -f pulsar.par pulsar.tim -theme dark -palette contrast

Themes are `dark`, `light` and `classic`; palettes are `ocean`, `contrast` and
`sunset`. The bottom menu's Light/Dark and palette buttons change them while
the viewer is open. Setup files accept `theme`, `palette` and `pointsize`
(a multiplier, default 1.3). Command-line theme and palette choices take
precedence. Explicit custom setup-file background/line colours are respected.

The modern layout adapts to window width and height. **Jumps** have
their own subarea below the toolbar: buttons wrap into columns, the window
grows when a longer list is selected, and `<` / `>` page through lists that
exceed the available space. Automatic growth is capped at 1000 plotting pixels;
manual resizing remains available. Jump buttons use contrasting text in both
selected and unselected states and toggle their fit flags.
Flag-selected jump labels show the group value alone, for example `UWL_PDFB4`
instead of `J -sys UWL_PDFB4`. Instrumental correction flags following the
convention `-<instrument>_<five-digit MJD>_jump 1`, such as
`JUMP -medusa_58925_jump 1 -2e-7 1`, appear as plain, non-clickable text with
a `fixed` annotation and the full selector name (`medusa_58925_jump`). There
is no button background or border. Their offsets remain applied, but `gplk` excludes
them from initial fits, refits and reload fits, and clicking the card cannot
enable fitting. Ordinary numeric selectors such as `-sys 1` remain selectable.
This is a `gplk` policy; the input parameter file is not rewritten. For other
TEMPO2 tools, explicitly end the correction line with `0` to keep it fixed:

    JUMP -medusa_58925_jump 1 -2e-7 0

Dedicated boolean selectors use the flag name: `JUMP -blahblah 1 8989.9898`
is displayed as `blahblah`, rather than `1`. Conventional grouping flags
(`-sys`, `-f`, `-g`, `-be`, `-backend`, `-group`) still display the selected
group value. Name extraction does not change which jumps are eligible for fitting.

The right-hand **Color by** panel offers a dropdown for frequency, a single color, setup-file
rules, observatory, and every flag ID found in the observations. Selecting a flag
also lists its values and color swatches in that panel, with its own paging.
Changing legend pages does not move jump buttons. Observations missing that
flag get a separate `(missing)` group. Long labels are abbreviated visually;
their stored identifiers remain unchanged. `Ctrl+I` retains terminal selection.

Automatic resize redraw requires a Giza X11 driver that sends resize
notifications through the motion-callback API. Older drivers, including the
tested Giza 1.5 build, can use **Shift+R** after resizing. Lowercase `r` keeps
its existing reload/refit meaning. Resize redraws preserve the fit, zoom and
overlays; the separate spectrum window supports the same redraw handling.

Drawing uses native Giza frame buffering and one point batch per selected
color group, with reusable storage. Flag coloring replaces frequency coloring
instead of drawing the same residuals again. The double-precision batch API
also avoids premature buffer flushes in some Giza float-point implementations.

Overlays use styles independent of the residual palette: binned points are
large gold (dark backgrounds) or burnt-orange (light backgrounds) diamonds,
with opaque error bars and contrasting outlines. Bins draw above other
overlays. Shapiro and FITWAVES curves use bold magenta dashed lines, and
selected points use outlined markers. A small caption identifies active
overlays. These styles also apply when using the classic residual palette.
Press `b` to compute bins; press `b` again to remove their overlay without
another width prompt. Pressing `b` once more computes a new set. The top-bar
**Clear overlays** button removes bins, model/solar/noise curves, selections,
labels and guide markers together, without refitting or changing subtraction
flags. `Ctrl+B` and `Shift+B` retain their older commands.

The **Shapiro** button toggles the stored solar Shapiro delay, in the current
residual units. `Ctrl+U` also toggles it and prompts for a vertical offset in
those units. This is a reference curve for a delay already included in the
timing model, not an extra correction applied to the residuals.

The **Solar** button adds geometric Sun–pulsar angular separation at the TOAs,
using the timing solution's observer, Sun and pulsar vectors. Its green curve
uses a separate 0–180 degree right axis, replacing the optional pulse-period
axis while active. This shows sampled solar approaches; it does not solve for
an exact conjunction between TOAs. These overlays use date/year/MJD axes.
They can also be enabled for file output:

    tempo2 -gr gplk -f pulsar.par pulsar.tim -us -solar -shapiro

For a fitted timing-noise model, use the Red, DM, Chrom, Total and Clean
buttons, press `W` to cycle the available components, or select a view directly:

    tempo2 -gr gplk -f pulsar.par pulsar.tim -noise red -us
    tempo2 -gr gplk -f pulsar.par pulsar.tim -noise clean -us

The reconstructed red, DM and chromatic signals are shown as timing delays at
each TOA's observing frequency, in the selected residual units. Total sums
these three components; Clean displays postfit residuals minus their sum.
These views do not alter observations, fit parameters or subtraction flags.
Use Residual to return to the ordinary postfit plot. Components without a
reconstruction are disabled, and a requested unavailable view reports why.
Individual components use their stored model errors. Total has no combined
error bars because cross-component covariances are unavailable. Clean shows
TOA measurement errors; these do not include the noise-model uncertainty.

These are fitted time-domain reconstructions, not newly generated stochastic
realizations or posterior draws. The top-bar **Noise overlay** button cycles
the available red/DM/chromatic/total reconstructions over the residuals, then
turns them off. Cyan hollow squares distinguish them from residual circles
and binned diamonds. Only achromatic red noise is joined by a dotted curve;
DM/chromatic/total values stay at their actual TOA frequencies, without joining
different observing bands. Individual components use stored uncertainties.
When residual mean removal is active, the overlaid reconstruction is centred
on the same plotted sample set. For example:

    tempo2 -gr gplk -f pulsar.par pulsar.tim -us -noise-overlay red

Click Spectrum or press `|` to inspect the current residual/noise series in
a separate window; `q` closes that window and returns to the residual viewer.
Use **Zoom** or `z`, then click two opposite corners to zoom in frequency and
power. `+` / `-` zoom around the view centre; **Reset**, `u`, or right click
restores the full spectrum. Esc cancels a selection. Resizing/redrawing keeps
the chosen zoom. Zoom changes the displayed limits, not the computed samples.

The spectrum toolbar offers three estimators:

 - **GLS** (`1`): weighted floating-mean Lomb–Scargle, normalized power 0–1.
 - **TEMPO2 WLS** (`2`): TEMPO2's native `TKspectrum` weighted least-squares
   estimator, with its native PSD normalization in yr³.
 - **Cholesky** (`3`): TEMPO2's native covariance construction and `calcSpectra`
   estimator, also in yr³. Select a covariance-model file using **Model** (`m`),
   or provide `-spectrum-model file`. An existing `-dcf file` is used by default.

Model and **Save data** (`s`) use a text field inside the spectrum window:
type the path and press Enter; Esc cancels. Save exports the full numerical
spectrum, including units and frequency columns in cycles/day and cycles/year,
and refuses to overwrite an existing file. Changing estimators resets the
display limits and does not refit the pulsar or reorder its observations.
Spectra can also be exported without opening the interactive viewer:

    tempo2 -gr gplk -f pulsar.par pulsar.tim -spectrum residual -theme light -grdev spectrum.pdf/pdf
    tempo2 -gr gplk -f pulsar.par pulsar.tim -spectrum red -spectrum-fmin 0.0001 -spectrum-fmax 0.05
    tempo2 -gr gplk -f pulsar.par pulsar.tim -spectrum residual -spectrum-method tempo2 -spectrum-output wls.tsv -grdev wls.pdf/pdf
    tempo2 -gr gplk -f pulsar.par pulsar.tim -spectrum residual -spectrum-method cholesky -spectrum-model covariance.model -spectrum-output chol.tsv -grdev chol.pdf/pdf

Spectrum sources are `residual`, `red`, `dm`, `chromatic`, `total` and `clean`.
All estimators use valid non-deleted TOAs within any active START/FINISH fit
range and their measurement errors, independently of the main plot's zoom.
GLS uses barycentric epochs; the native methods use site epochs, matching
TEMPO2's spectral tools. The default weighted, floating-mean GLS calculation's
standard power is the fractional reduction in weighted squared residuals
relative to a constant model, between 0 and 1. The dashed curve is the weighted
sampling window, useful for spotting cadence aliases. The strongest sampled
peak is labelled with its period in days. This diagnostic is not a calibrated
power spectral density or a detection-significance calculation.

For GLS, default bounds run from one cycle per data span to a cadence-based upper
bound; irregular sampling has no universal Nyquist limit. Override them with
`-spectrum-fmin` and `-spectrum-fmax` in cycles/day. The frequency grid samples
at least five times per inverse data span where possible, with 256–8192 bins;
narrow the bounds to inspect fine features over long baselines. The plot's
frequency axis is in cycles/year. The method and normalization follow
[Zechmeister & Kürster (2009)](https://doi.org/10.1051/0004-6361:200811296) and
the [Lomb–Scargle normalization documentation](https://docs.astropy.org/en/stable/timeseries/lombscargle.html#periodogram-normalizations).
Native methods keep TEMPO2's Fourier grid; frequency bounds filter those bins
without oversampling them. Native covariance-model syntax and validation apply
to Cholesky spectra. These diagnostics do not supply false-alarm probabilities
or posterior uncertainty bands.

5.2 Building your own plugin
----------------------------
The easiest way to compile your own plugins is:

    g++ {$CFLAGS} {$LDFLAGS} -fPIC -shared -o {$TEMPO2}/plugins/{$PLG_NAME}_{$LOGIN_ARCH}_plug.t2 {$SRCLIST}

where:

 - `{$PLG_NAME}` is the name of your plugin
 - `{$SRCLIST}` is your plugin's source code.
 - `{$LOGIN_ARCH}` is the result of `` `uname` `` (usualy Linux).
 - `{$CFLAGS}` are the compiler flags your plugin needs... remeber to add a -I option to point to the location of tempo2.h
 - `{$LDFLAGS}` are any linking options you need, e.g. pgplot, etc.
 - `{$TEMPO2}` is the tempo2 runtime dir

For example, to compile a basic plugin called '`foo`' on linux, you might do

    g++ -I/usr/src/tempo2 -fPIC -shared -o $TEMPO2/plugins/foo_{$LOGIN_ARCH}_plug.t2 foo_plug.C


5.3 Adding a new plugin to the default build list
-------------------------------------------------
If your plugin has dependances that are already covered by the lists above,
just add the name to the appropriate list, and name your plugin source file as:

name_plug.C

6. Changes from the old Make system.
------------------------------------
At the start of 2010, tempo2 moved over to an autotools based make system,
replacing the old hand written makefiles. This may confuse some people!

Important notes:

 - Tempo2 plugins now have a .t2 extention, rather than the old .so
    This is to ensure reduce confusion on MacOSx and to allow the old
    make system and the new make system to co-exist for a while.
 - Any 3rd party plugins will still work as before. Indeed, to update
    a plugin, just change the .so extention to a .t2 extention.
    e.g. mv general_Linux_plug.so general_Linux_plug.t2

7. Installation Troubleshooting
-------------------------------

7.1 Can't find PGPLOT
---------------------
Download pgplot from:
http://www.astro.caltech.edu/~tjp/pgplot/

Or use PSRSOFT to manage the installation.
http://www.pulsarastronomy.net/wiki/Software/PSRSoft

If you have pgplot installed, but it is not detected by the configure script, check:

 - You have got at least libpgplot.a and libcpgplot.a in your LDFLAGS
 - Check you have `$PGPLOT_DIR` pointing to the folder with grfont.dat and rgb.txt
 - Check that you have `$F77` set to the same compiler that compiled PGPLOT
    (e.g. setenv F77 gfortran, if you used gfortran for PGPLOT)


7.2 Incompatible C and Fortran compilers
----------------------------------------
Check that you are using the same build of gcc and gfortran (or whatever compiler you are using).

Note that on MacOSX there is often an issue where the default compiler is incompatible with gfortran.
The gfortran compatible version is often called gcc-6 and gxx-6 or similar. Use this with:

    export CC=gcc-6
    export CXX=g++-6

and reconfigure.


8. Bugs and feature requests
-----------------------------
Please submit bug reports here: https://bitbucket.org/psrsoft/tempo2/issues/new

Note that it is very helpful if you can upload a small example demonstrating the bug!
