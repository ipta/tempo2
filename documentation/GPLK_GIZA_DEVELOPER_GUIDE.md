# Maintaining the native Giza gplk plugin

Start here when changing `tempo2 -gr gplk`. The plugin calls the native Giza C
API and links `libgiza`. The compatibility-named palette constant
`GIZA_COLOUR_PALETTE_PGPLOT` only selects colour values; it does not introduce
a dependency on `libpgplot` or `libcpgplot`.

## Source map

Paths below are relative to the repository root.

| File | Responsibility / useful entry points |
| --- | --- |
| `plugin/gplk_plug.C` | Plugin entry point, argument/setup parsing, fitting, plot arrays, event loop and legacy commands. Start with `graphicalInterface`, `doPlot`, `setPlot`, `setLabel`, `drawMenu3_2` and `checkMenu3`. |
| `plugin/gplk_giza.h` | Device setup, themes, palettes, residual/overlay drawing, viewport conversion and cursor input. |
| `plugin/gplk_ui.h` | Bottom jump buttons, right-hand color panel/dropdown, dynamic color-group cache and literal identifiers. |
| `plugin/gplk_jumps.h` | Concise selector labels, fixed instrumental-correction convention, fit-policy enforcement and guarded toggling. |
| `plugin/gplk_noise.h` | Read-only access to reconstructed noise signals, availability checks and uncertainties, all in seconds. |
| `plugin/gplk_overlays.h` | Explicit overlay context, solar-angle geometry, current-unit Shapiro/noise overlays and noise-overlay cycling. |
| `plugin/gplk_spectrum.h` | Spectrum sample/result types, numerical calculation (`gplkComputeSpectrum`) and device/UI integration (`gplkShowSpectrum`). |
| `plugin/gplk_spectrum_native.h` | Native TEMPO2 WLS/Cholesky adapters, fit-span filtering, stable row/uncertainty/index ordering. |
| `tests/testGplkGiza.C` | Headless renderer integration, noise units/availability, and spectral regression tests. |
| `autoconf.boot/giza.m4` | Native API link probe and `GIZA_CFLAGS` / `GIZA_LIBS`. |
| `bootstrap`, `plugin/plugin_lists/giza.plugins`, `plugin/Makefile.am` | Generate and select the Giza plugin target independently of PGPLOT. |
| `README.md` | User-facing build instructions, options and scientific interpretation. |

The plugin build selects `plugin/gplk_plug.C`. Files in `plugin/plkfast/` are
not selected by this build; do not assume their similarly named sources are
the active implementation.
The original PGPLOT implementation remains `plugin/plk_plug.C`, built from
`plugin_lists/pgplot.plugins`; native `gplk` is selected by `giza.plugins`.
Native setup files are `T2runtime/plugin_data/gplk_setup*.dat` and are installed
by `gplk-install`. Keep Giza names and settings separate from the original plugin.

The helper headers deliberately use internal-linkage functions and state in
the existing C++03-compatible, single-source plugin. They are also included by
the test executable, which gets independent state. If splitting the plugin
into multiple translation units, first move shared configuration into one
implementation file or an explicitly passed context; otherwise each source
file will silently get its own theme, palette and spectrum bounds.

## Execution and data flow

`graphicalInterface` resets per-invocation display settings, parses options,
loads/preprocesses the pulsar and calls `callFit`. A standalone `-spectrum`
request then goes directly to `gplkShowSpectrum`; otherwise execution enters
the `doPlot` loop.

`doPlot` reads setup options, opens a device and repeatedly:

1. Builds plotted `x`, `y`, lower/upper uncertainty and observation-index arrays.
   `setPlot` selects a quantity and `setLabel` describes its units.
2. Applies existing mean/phase/display conversions, sorts all parallel arrays,
   and establishes the viewport and world-coordinate limits.
3. Buffers a frame: panels/grid/axes, residuals and their subdued error bars,
   model curves and selections, then bins and the overlay caption.
4. Ends buffering, reads a cursor/key event and updates the next frame's state.

`id[i]` maps each plotted sample back to `psr[0].obsn[id[i]]`. Keep it aligned
with every reordered data/error array; deletion and identification depend on
this mapping. Display helpers receive borrowed arrays and must not modify
their contents or the pulsar's fitted values.

## Themes, palettes and drawing state

`gplkTheme` chooses dark, light or classic appearance. `gplkPalette` chooses
residual-series colours. Explicit command-line choices take precedence over
setup-file choices. Without an explicit theme, `-publish` / `-image` choose
light and normal viewing chooses dark. `doPlot` also replaces the old default
runtime setup's font/series styling for modern themes, so an older installed
`gplk_setup.dat` does not force the previous appearance. Explicit custom setup
files retain their supplied styling.

The colour table has two separate purposes:

| Indices | Meaning |
| --- | --- |
| 0, 1 | Canvas background and foreground. |
| 2–12 | Residual series; the selected six-colour palette repeats. |
| 13–15 | Deleted-point / muted legacy colours. |
| 256–260 | `GPLK_PANEL`, `GPLK_GRID`, `GPLK_MUTED`, `GPLK_ACCENT`, `GPLK_ERROR`. |
| 261–263 | `GPLK_BINS`, `GPLK_MODEL`, `GPLK_SELECTED`; overlay roles independent of the series palette. |
| 264–265 | `GPLK_SOLAR` and `GPLK_NOISE_OVERLAY`, separate from series colours. |

Keep these reserved indices out of user series assignments. Check the
supported Giza colour-index range before reserving more entries.

Giza is stateful. Drawing helpers generally bracket attribute changes with
`giza_save()` / `giza_restore()`. These restore drawing attributes, **not**
viewport/window geometry, device selection or colour-table definitions.
Restore those explicitly if a new helper changes them. Updating reserved
scratch/role colours is intentional; overwriting a residual-series entry is
not. UI helpers assume the caller has already selected the appropriate
viewport and world-coordinate window.

Residual error bars go through `gplkDrawErrorBars` (alpha 0.42, width 0.65 in
modern themes). Bins go through `gplkDrawBinnedOverlay`: opaque error bars,
background halos, and outlined diamonds. Its nominal six-pixel marker radius
is converted to world-coordinate offsets, keeping bins legible across zoom
levels. Model curves use `gplkDrawOverlayLine` and selections use
`gplkDrawHighlightedOverlay`. `gplkPrepareOverlayColours` owns their colours;
`gplkDrawOverlayLegend` owns their caption. Bins draw last to remain visible.

## Coordinates and input

Native `giza_band` returns **zero on success**. Use `gplkReadBand`, which
returns a C++ boolean, preserves coordinates on failure and supplies a quit
key. Callers finish buffered drawing before waiting for input. Do not treat a native
nonzero return as a successful event.

Cursor positions are in the current window's **world coordinates**.
`gplkGetViewportPixels` derives bottom-left-origin physical coordinates from
the normalized viewport and paper dimensions. This works around reflected
vertical physical-viewport queries in the Giza version used during porting.
Use the same conversion for new hit testing and drawing; do not independently
invert the y axis.

The modern bottom toolbar is drawn by `gplkDrawFooter` in `gplk_ui.h`.
`gplkHandleUi` dispatches the right-hand color panel/dropdown and bottom jump subarea; `checkMenu3` dispatches
the original toolbar actions. The toolbar has two 24-pixel rows with ten
columns; the independent jump subarea wraps buttons below it:

| Row | Columns 0–9 |
| --- | --- |
| 0 | Quit, Clear, Measure, Help, Red, DM, Chrom, Total, Clean, Residual |
| 1 | List, Print, Highlight, theme toggle, palette cycle, Spectrum |

Change drawing and hit testing together. `GplkLayout` in `gplk_giza.h` computes
the top, left sidebar, plot, right color panel and bottom jump rectangles from surface pixels;
both renderers and hit tests use these same rectangles. The toolbar's physical
height stays fixed while jumps wrap according to width. `gplkFooterCell` owns
both button geometry and hit boxes, and the jumps header owns paging buttons.
`GplkUiState::jumps` contains only jump/FD-jump buttons; `legend` contains only
color-group entries. Jump text uses foreground color 1 on inactive buttons and
background color 0 on active accent buttons, independent of series colors.
`gplkDrawColourPanel` draws a separate right-side panel with colored swatches
and foreground legend text. `gplkColourSelector` owns the dropdown hit box;
the dropdown opens down over the legend. `page` and `legendPage` are independent,
so changing grouping or legend pages cannot move or toggle jump buttons.
Jump labels omit the `J` and flag-selector prefixes for ordinary flag groups;
MJD/FREQ selectors retain both bounds, and FD/DM jumps retain their type/index.
`gplkDescribeJump` is the shared source of labels and correction classification.
Only boolean-valued selectors named `-<instrument>_<five-digit MJD>_jump 1`
are classified as fixed corrections. Do not broaden that to every non-`sys`
flag or every selector whose value is `1`; those can be ordinary fitted groups.
Corrections render as text only, with a small `fixed` annotation and the full
selector name minus its leading dash. Do not draw a fill, outline or button
for these entries, and ignore their click area. Their offsets and TOA memberships remain
unchanged. `gplkToggleJumpFit` guards both modern and legacy mouse dispatch.
Label extraction is separate from fit policy: dedicated boolean selectors
such as `-blahblah 1 <offset>` use `blahblah` as their name, whether or not the
selector matches the fixed-correction convention. Conventional grouping flags
(`-sys`, `-f`, `-g`, `-be`, `-backend`, `-group`) retain their matched group value,
including numeric group names such as `1`.

TEMPO2 treats a trailing `1` as a fit request, including for these corrections.
`gplkKeepCorrectionsFixed` therefore runs before `callFit` and `reFit`, covering
startup, reload and replay paths; it clears only matching `fitJump` entries.
Never zero `jumpVal`, change `jumpSAT`, remove a selector, or clear a correction's
stored uncertainty to implement this policy. `formResiduals` applies offsets
independently of fit flags. The native parser and the on-disk parameter file are
unchanged; other TEMPO2 tools need an explicit trailing `0` for fixed jumps.
Growth is requested once on initial display or a new color selection, capped
at 1000 plotting pixels. Never repeatedly impose paper size during a user's
resize drag. Classic menus and explicitly configured non-menu viewports retain
their existing geometry.

`GplkColourGroups::rebuild` discovers flags and assigns observation IDs to dynamic
groups, including missing values. Rebuild after replacing observations; display
positions are separately mapped through the current `id[]` array. `groupPositions`
partitions points once per frame so drawing costs O(TOAs + groups), rather than
rescanning every observation for each flag value. Each group draws its bars and
points once; no frequency pass precedes categorical coloring. Scratch arrays
and conversion vectors reuse their allocations. Setup rules retain their legacy
overplot semantics. Color indices cycle through 2..12 for larger category sets.
Flag values are labels, not Giza math strings: `gplkLiteralText` preserves literal
underscores, braces, carets and backslashes in identifiers. Ellipsis changes
only the displayed label.
`gplkTextSize` queries native pixel text lengths and maps distances through the
current viewport/window; Giza's world-unit text-length conversion is incorrect
in some builds and can place the rest of a literal identifier off-screen.

`gplkReadEvent` is the idle main/spectrum input path. A motion callback opts into
resize notifications on newer Giza X11 drivers. These drivers return synthetic
`r`, so compare surface dimensions before/after input and return a redraw-only
event (key zero) when size changed. Never let this event reach lowercase `r`'s
reload/refit command. Shift+R explicitly redraws on older Giza drivers; tested
Giza 1.5 lacks automatic notification. Synchronize old-driver surface geometry
with `giza_change_page` on initial display/resize. File devices also need an
initial synchronization before the first drawn frame: otherwise requested
paper dimensions can differ from the actual surface and clip the jump row.
Do not advance file pages again during ordinary frame redraws.
Ordinary frames clear their background without recreating the device surface.
Rubber-band operations retain `gplkReadBand`, without opting into idle resize
events mid-selection. Their completed input returns to the normal redraw loop.

Each full frame is bracketed by native `giza_begin_buffer` / `giza_end_buffer`.
Use `gplkDrawPoints` for float arrays: it reuses double conversion storage and
calls native `giza_points`. The Giza source tested here ends buffering
unconditionally inside `giza_points_float`, causing partial-frame presentation;
the double entry point preserves it. `nativeRenderer` tests this property.
Do not add an extra flush just before input. Giza's public API does not provide
retained UI layers or a dirty-region cache: this implementation still redraws
a complete frame on state changes, and does no drawing on idle pointer motion.

Noise mode arrays appear in drawing/dispatch; keep their ordering aligned.
Keyboard dispatch lives in `doPlot`:
`W` cycles available noise modes, `|` opens a spectrum, and the older letters
already have many assignments. Search existing key handling before adding one.

The top action row is drawn by `drawMenu` and dispatched by `checkMenu`:
columns 4–5 clear overlays, 6 toggles Solar, 7 toggles Shapiro, and 8–9 cycle
noise overlays. `GplkOverlayControls` is owned by `graphicalInterface`, passed
to `doPlot`/menus, and holds these new layers explicitly. A clear request is
consumed at the start of the next frame, where the viewer owns all legacy
overlay flags. It does not change `TNsubtract*`, residuals or fit parameters;
it resets only internal selection flags and visual/cached state.

`gplkSolarSeparation` forms observer-to-Sun as `sun_ssb - earth_ssb -
observatory_earth`, with the same `STL_FBAT` spacecraft convention as
`shapiro_delay.C`. Angular separation uses normalized dot/cross products and
atan2, returned in degrees. The plotted solar curve is mapped to a distinct
0–180 degree secondary axis; it must never be added to residuals or interpreted
in seconds. Invalid/unavailable geometry is rejected. This is geometric
separation at sampled epochs, not an apparent-position or conjunction solver.

Shapiro/noise delays use `gplkOverlayDelay` to match `setPlot`'s current units
and phase convention, and borrow the current sorted x/id arrays instead of
assuming MJD-minus-epoch or microseconds. Joining DM/chromatic samples across
observing frequencies would imply a misleading continuous curve; those
reconstructions use separate markers/errors, while red noise can be joined.

`gplkShowSpectrum(..., hasParent=true)` saves the parent device ID, opens its
own device, waits for close input, closes only its own device and reselects
the parent. A standalone call uses `hasParent=false`. Do not replace this
close with `giza_close_devices()`, which would also close the main viewer.

## Noise quantities and units

The added modes form a contiguous range used by dispatch, scaling and menus:

| Mode | ID | Value returned by `gplkNoiseValue` |
| --- | --- | --- |
| Red | 20 | `TNRedSignal` |
| DM | 21 | `TNDMSignal` |
| Chromatic | 22 | `TNChromSignal` |
| Total | 23 | Sum of the three configured components |
| Clean | 24 | Raw postfit `residual` minus that sum |

These reconstructed signals are already timing delays in **seconds at the
TOA's observing frequency**. Do not apply another DM/chromatic frequency
factor. The older modes 16/17/19 have different legacy representations and
are not interchangeable with these modes.

Always check model configuration before reading component fields: the core
does not initialize absent components universally. `gplkNoiseAvailable` also
requires a non-deleted sample with nonzero signal or uncertainty. Zero or
unavailable reconstruction is not proof that the pulsar has no timing noise.

`setPlot` divides the new modes' seconds by `unitFlag` for display and applies
the existing phase conversion when requested. The `doPlot` uncertainty path
expects microseconds before its common display-unit conversion, so it
multiplies `gplkNoiseUncertainty`'s seconds by `1e6`. Do not scale twice.

Individual components use the stored model standard errors. Total returns no
combined uncertainty because component cross-covariances are not available.
Clean uses TOA measurement errors and does not propagate noise-model errors.
Neither operation alters fitted parameters, observations or `TNsubtract*`
flags. Group-noise fields exist in the structure but this core does not
populate a reconstruction for the new viewer; a new component needs an actual
producer, initialization and validity checks before it can be displayed.

## Spectrum calculation

`gplkComputeSpectrum` (GLS) is independent of graphics and receives borrowed sample
values through a const vector. Sample times are days, values/sigma have the
same units (seconds in the plugin), and input/output frequencies are
**cycles/day**. Only plotting converts frequency to cycles/year and log10
coordinates. Peak periods are reported in days.

The calculation rejects invalid samples/nonpositive errors and requires at
least five samples spanning multiple epochs. It fits an intercept plus sine
and cosine at each trial frequency using inverse-variance weights. Common
rescaling of all sigma values cancels. Weighted centering of sine/cosine
implements the floating intercept, and a determinant check skips singular
fits. Output standard power is `(chi2_constant - chi2_fit) / chi2_constant`.
The sampling-window curve is `|sum(w * exp(i * phase))|^2 / sum(w)^2`.

Default bounds are one cycle per span and a median-positive-cadence upper
bound; the latter is not a universal Nyquist frequency. The uniform trial grid
targets five samples per inverse span, bounded to 256–8192 bins. At the cap,
a narrower requested frequency range is needed for finer peak localization.
Complexity is O(number of TOAs × number of trial frequencies), with no FFT or
external numerical dependency. Powers are diagnostic, not a calibrated PSD
or false-alarm significance. The method citation is in `gplk_spectrum.h`.

`gplk_spectrum_native.h` provides two separate native PSD adapters. TEMPO2 WLS
calls `TKspectrum` with specType=4, specOut=1, no interpolation/prewhitening,
ofac=hifac=1. Cholesky calls `getCholeskyMatrix` and `calcSpectra` with nfit=-1.
Do not replace either with the GLS implementation or rescale their PSDs to 0–1.
Native output is yr^3; x remains cycles/day. These paths use site epochs minus
52601, matching `cholSpectra`; GLS uses barycentric epochs. All paths apply the
same deleted/START/FINISH selection, including BAT-based bounds (fitFlag=2).

Sort native time/value/error/observation-ID tuples together before invoking
TKspectrum: its internal sorter only permutes x/y. Native Cholesky receives
a metadata copy because covariance construction writes detL. Observations
are borrowed read-only, and any extra covariance is projected into the sorted,
filtered row order. The native covariance-model parser retains TEMPO2's syntax
and error handling. A missing/unreadable model is rejected before entry.
Dense Cholesky storage/cost follows the native implementation and can be
substantial for large datasets. Tests check a known sinusoid's absolute PSD,
nonzero covariance, native WLS agreement, filtering, and fit preservation.

`GplkSpectrumView` owns display bounds only: log10(cycles/year) and either GLS
power or log10(PSD / yr^3). Rectangle selection clips to the current view;
plus/minus scale around its centre, and reset restores saved full bounds.
Resize redraws reuse those bounds and the computed result. Method changes
compute a trial result first, then atomically replace the curve, normalization
and full bounds. Zoom never recomputes the estimator or changes the TOA set.
The sampling-window curve is shown only for GLS, never mixed onto PSD axes.

Spectrum controls use `gplkSpectrumButton` for drawing and hit testing. The
in-window Model/Save text prompt redraws only its input area, then restores the
plot geometry. Numerical export writes the full computed spectrum with method,
series and units; exclusive creation prevents silent replacement of a previous
result. File-device plots omit interactive buttons. Parent device handling
remains in `gplkShowSpectrum`; never close all devices to leave this child.

## Common changes

**Change bin colour, marker size, outlines or drawing priority:** edit
`gplkPrepareOverlayColours` / `gplkDrawBinnedOverlay` in `gplk_giza.h`; priority
is at the end of the buffered frame in `doPlot`. Keep residual and overlay
styles distinguishable by shape/width as well as hue. This does not require
changing bin statistics.

**Add a palette:** extend `GplkPalette`, the parser and the
`colours[palette][light_or_dark][six_series]` table. Also update the toolbar
names, its cycle modulus, CLI help and README. Check dark/light and custom
setup-file precedence, not just a single export.

**Add an overlay:** calculate its values separately, draw in current world
coordinates with saved/restored attributes, reserve a semantic style if
needed, then add the active caption. Place it before bins unless it should
intentionally cover them. Preserve plotting-array/index alignment.

**Add a noise component/view:** provide initialized, meaningful core fields;
extend mode parsing, names, configuration/availability, values and uncertainty
policy. Search all uses of `GPLK_NOISE_CLEAN` and the existing mode arrays:
range checks, keyboard wraparound, `setPlot`, `setLabel`, error scaling,
toolbar rendering/dispatch, CLI help and spectra all need consistent handling.
Add numerical tests for units, missing data and subtraction before UI testing.

**Change spectrum methods or normalization:** keep calculation separate from
device handling, update the formula/axis labels and README together, and
compare against an independent weighted linear least-squares calculation.
Retain coverage for uneven sampling, heteroscedastic errors, offsets, common
error scaling and degenerate inputs.

## Build and verification

From a configured source tree with Giza:

```sh
./bootstrap
./configure --with-giza    # or --with-giza=/installation/prefix
make
make -C plugin gplk
make -C tests check TESTS=testGplkGiza check_PROGRAMS=testGplkGiza
```

Set `TEMPO2` to an installed runtime directory before configuring/running;
use `--prefix` and `--with-tempo2-plug-dir` when selecting install paths.
For a local uninstalled plugin, `TEMPO2_PLUG_PATH="$PWD/plugin/.libs"` selects
it when invoking the matching built `./tempo2`. Check the reported plugin path
so tests do not accidentally load an older installed build.

The tests exercise headless SVG output, coordinates, noise data access and
the spectrum's independent numerical reference. They do not replace a desktop
check of mouse input, window resizing or controls. For rendering changes,
inspect dark/light outputs with dense residuals, bins, selections and a model
curve together. For scientific changes, add focused numerical tests.

To verify native-only graphics linkage on ELF platforms:

```sh
readelf -d plugin/.libs/gplk_*_plug.t2
nm -D plugin/.libs/gplk_*_plug.t2
```

Expect a `libgiza` dependency and native `giza_*` references; investigate any
`libpgplot` / `libcpgplot` dependencies or `cpg*` symbols. Include a new helper
header in `plugin/Makefile.am`'s `EXTRA_DIST` so source distributions retain it.

## Existing limitations to keep in mind

The old `gplk_plug.C` remains a large function-oriented implementation with
manual allocations, fixed-size arrays and duplicated menu tables. This guide
documents the new boundaries; it does not mean the whole legacy plugin has
been refactored. `viewModels` is explicitly disabled in this source;
Shapiro/FITWAVES drawing is separate and remains available.

Binning retains its legacy inverse-error weighting and cached overlay
behaviour; it is not the inverse-variance weighting used by the spectrum.
Recompute bins after changing the data being compared. Changing bin statistics
or cache invalidation should be a separate scientifically validated change.
The main-loop `b` key clears an existing bin count/redraw flag/cache position,
or computes bins if none exist. Automatic overlay redraw requires a positive
bin count. Clear overlays resets the same cache along with other visual layers.

Spectra use valid non-deleted raw per-TOA values within active START/FINISH fit bounds, not the current zoom,
highlight selection, averaged display arrays or the ordinary plot's
`residualtn` substitution. Clean is its own explicit subtraction view.
Standalone spectra also bypass `doPlot`'s setup-file parser; use the CLI for
their theme/palette settings. New theme/noise/menu controls are not all
serialized by legacy record/replay; add matching reader/writer support if a
feature needs reproducible visual replay.
