Koch-prefractal waveguides
Research data accompanying the manuscript

The cell has period P = 2 pi, plate separation b = pi and c = 1.
Frequencies and Bloch wavenumbers are reduced. Columns ending in _pct
are percentages. qI and qEM in fields/ are fractions of the whole-cell
integrals unless a column says otherwise. CSV files use UTF-8, period
decimal separators, and a header. The string nan denotes an unavailable
crossing or estimator. Retain it as missing.

geometry/unit_side_vertices.csv
Closed polygon contours for the six base shapes at orders 0 to 3
and for the four regular polygons at order 4.
Multiply a unit contour by L and translate it to the cell center
(pi, pi/2). The FK0 contour is shared by the two growth directions.
geometry/fk3_fk4_surfaces/ contains the collocation points, panel
lengths and normals used for the 24 inward FK3/FK4 comparisons,
including the two rectangle orientations.

gaps/hull_comparisons.csv and gaps/inward_edges.csv
Numbers behind the comparison of conducting Koch prefractals
with their hulls and same-area twins. The three base-polygon filling
levels in inward_edges are distinct from each row's actual filling.
edge_shift and sigma are reduced frequencies. The relative differences
in hull_comparisons are 100 times (prefractal/reference - 1).
gaps/equal_filling.csv contains the actual equal-area series and the
widest relative gap for each feasible shape, order and orientation.

gaps/smooth_curves.csv, saving_methods.csv, paired_curves.csv,
paired_factors.csv and paired_runs.csv
The dense curves give the widest full gap below omega_r = 2.6 as a
function of smooth filling. saving_methods contains the 36 inward target
cases. The twelve individual method columns are interpolator, window,
subsampling and calibration variants. shifted_low and shifted_high change
the target by 0.2 percentage points in either direction using the central
quadratic estimator. The paired data use five meshes: E standard,
G inclusion refinement, P plate and periodic-cut refinement, R spectral
refinement, and J all three. Each mesh has one inward triangular FK3
case and eleven smooth triangles. Stored paired roots in bands/paired/
have two whitespace-separated columns, K_r and omega_r; each case name
matches paired_runs.csv. The full-precision band roots are used by the
checks below. The dense smooth curve is rounded to 0.001 gap percentage
point in its archived table.

bands/*.csv
The six individual root lists cover the triangular controls of the
fixed-side and equal-filling comparisons. square_parity contains
three square frequency windows with two branches and 51 K_r samples per
branch. seed_window is the order and seed band, parity is the measured
reflection in y. A nonzero subspace_flag denotes treatment of close
roots (2 marks the second-singular-value search). The triangle_modes
table contains 56 field-matched narrow and dispersive bands, their
dense-trace widths, maximum slopes and the geometric attenuation T.
triangle_joint_traces and triangle_widths give the four inward FK3
and hull pairs at the three boundary refinements. Widths are in reduced
frequency units, including where the paper displays them in thousandths.

fields/
recess_by_mode.csv gives 214 zone-center cases containing recesses,
with volume, Green, midpoint, and absolute half-spread at their finest
available boundary resolution, as well as the two standard-mesh
estimates and the standard-mesh intensity fraction. The 21 modes with
best_mesh equal to 4 have controls at factors 1, 2 and 4. The other 193
use factor 1. recess_refinement.csv records the underlying volume
integrals on h/2, h/4 and h/8 and the Green-identity estimates. The
reference spacing is h = pi/120. pockets.csv gives the regional
weights where evaluated. Source tags use interna and externa for inward
and outward growth, igualn for rescaling to equal filling, and ILtriangulo
and ILpentagono for the 2.19% and 8.72% base-polygon levels.

triangle_maps.npz carries a subset of six 480 by 240 intensity maps,
their reduced frequencies, grid-center coordinates, and polygon vertices.
Keys twin, FK3 and hull identify the same-area smooth triangle, inward
triangle and its convex hull. The two row numbers are the zone-center
mode indices. Array indices are [y, x]. The intensities come from the
IEM field reconstruction at the recorded frequency and source density.
boundary_sources contains those six complex density vectors and the
boundary nodes and quadrature weights for each geometry. In mode_*.dat
the final two columns are Re(f) and Im(f). The fields displayed and the
boundary-refined energy estimates are identified separately.

manuscript_maps.npz contains the intensity arrays and exact polygons
used for the four multi-panel field figures. manuscript_maps.csv
records the frequency and final plotted weights for each field key.
manuscript_map_index.csv maps those keys to the mode and role in each
field comparison. Each intensity array has a key ending in _I.
intensity_correlations.json gives the intensity correlations in the
region common to the three geometries outside the hull.

validation/
fem_pares_TE.dat gives total widths and relative shape differences with
the IEM and two independent finite-element meshes. fem_vs_iem_TE.dat
lists the physical gap counts for the four benchmark geometries in the
full sweep and up to the analysis limit omega_r = 2.6 (Table 2).
timing.dat contains the elapsed times in seconds and peak memory in
kilobytes for the comparison reported in Table S2. Its variant tags
give the panel count N5 or FEM element size h and the thread count.
convergencia_resumen_TE.dat and
standard_and_refined_gap_edges.dat record contour and sweep controls.
convergencia.dat lists the inward refinement settings; inward_gap_edges.dat
contains the matching gap edges. Original column descriptions are in the
headers. In these source tables, prohibidas means band gaps, paredes means
plates and cuts, resolucion means spectral resolution, and nu is the number
of panels per segment.

validation/dispersion/
The eight point tables give K_r and omega_r for the triangle FK0,
inward triangle FK3, outward vertical rectangle FK3, and hexagon FK0.
The input spectra are in validation/dispersion/source_spectra/.
The FEM spectra use h = 0.04. To export the point tables from the
spectra, run python3 -B figures/prepare_band_panels.py OUTPUT.

validation/case_b_mesh/
The three boundary_m*.dat files are the stored surfaces for the
triangular hull in case B of Table S17. The first line gives the number
of contours. The second gives the panel counts for the two plates,
the two periodic cuts, and the hull, in that order. The remaining rows
give the panel midpoint, tangent, second derivative, and length.
The maximum hull panel length is 1/(11.5*m). A side of length L
uses ceil(11.5*m*L) equal panels. For m = 1, 2 and 4 the
files contain 54, 105 and 207 hull panels, respectively.

bands/inventory/
classified_bands.csv contains the fixed-side and equal-filling sets at the stated
frequency window. The smooth triangle FK0 at L = 1 is shared by the
two ensembles. Retain one occurrence for the 479-branch distribution.
The fixed-side set includes 29 branches in five FK4 structures.
endpoint identifies the branches that use a parity trace, a fitted
endpoint, or the resolved part of a trace. prohibidas_L1n2_TE.dat gives
the sweep gap edges used for shading the fixed-side plot.

bands/fk3_fk4/
The sweep files and 11-point traces reproduce gap edges for the 24
inward FK3, FK4 and hull structures at 2.19%. gap_edges.csv gives
the intervals used to compute total gap widths in that comparison.
fk4_width_controls.csv records the 113 seeded candidates, and
square_branches.csv gives all the parity-resolved square cases.
The physical flag excludes the interior resonance of the outward
FK4 hexagon. The odd FK4 square branch has a fit at K_r = 0 because
its minimum could not be separated from the even branch at that point.
The inward fixed-side triangular FK4 band 6 overlaps another root at
K_r = 0. Its width column gives the span of the ten resolved samples.
The archived widths and singular-value diagnostics remain in the
trace headers and summary files. Surface paths in the trace headers
are relative to the original IEM_resultados archive.

bands/inward_fk3_hull/
paired_widths.csv includes the 66 pairs reported in the article. The
resolved column identifies them among the candidate pairs. Trace files
use 11 wavenumbers for m1 and the two zone endpoints for the joint
boundary refinement m2.
The geometry (triangulo, cuadrado or pentagono) at the start of
each filename identifies the base filling level. Physical shape is
the next part of the tag. b91 labels one branch recovered by a local
root scan. resolved is true when both widths reach 1e-5 at both
refinements.

The folders gaps/shape_and_iteration_tables and bands/diagnostics
provide supporting values for the extended shape comparisons,
equal-aperture test, last resolved iterations and side thresholds.

Reproducing the reported comparisons
Run python3 reproduce_checks.py from the dataset root (NumPy required).
It reads the stored spectra through the shared gap detector, checks
all 60 paired largest gaps, fits their smooth reference curves, and
verifies the reported branch and energy results. For quick-look plots of the archived field and band subsets, run
python3 fields/plot_maps.py triangle_field_maps.pdf or
python3 bands/plot_triangle_bands.py triangle_bands.pdf.
The commands below reproduce the included publication figures.
Run python3 reproduce_extended_bands.py to check the branch
classification, FK4 widths, traced gap edges and FK3/hull bandwidth ratios.

Figures and plotting
figures/ contains the 13 main figures and five supplementary figures
as PDF files. From the dataset root, install the Python dependencies and
render the full set into a separate directory.

    python3 -m pip install -r requirements.txt
    python3 reproduce_figures.py ../reproduced_figures

The scripts use the bundled spectra, field arrays, geometries, and
derived tables. The commands below also generate individual figures or
groups.

Fig. 1  figures/waveguide_geometry.pdf
  Command  python3 figures/plot_geometry.py OUTPUT
  The unit-cell diagram and inclusion outlines use the stored
  vertices in geometry/unit_side_vertices.csv.

Fig. 2  figures/finite_element_band_validation.pdf
  Command  python3 figures/finite_element_band_validation.py --pdf OUTPUT/finite_element_band_validation.pdf
  The comparison uses the IEM and FEM point tables in validation/dispersion/.
  Red shading shows the FEM gaps up to omega_r = 2.6. Green markers
  identify complete IEM columns with sampled width W < 0.01.

Figs. 3–5, 13 and S1
  Command  python3 figures/plot_band_figures.py OUTPUT
  Fig. 3 uses bands/inventory/prohibidas_L1n2_TE.dat and
  bands/inventory/classified_bands.csv.
  Figs. 4 and 5 use the fixed-side triangular traces in bands/,
  bands/twin_outward_FK3_L1.dat and geometry/unit_side_vertices.csv.
  Fig. 13 uses bands/inward_fk3_hull/paired_widths.csv and plots the
  66 retained pairs. Fig. S1 uses bands/square_parity.csv and the two
  outward-square FK4 traces in bands/fk3_fk4/traces/.

Figs. 6–8 and S2
  Command  python3 figures/plot_control_figures.py OUTPUT
  Fig. 6 uses gaps/shape_and_iteration_tables/curvas_igual_n_TE.dat.
  Fig. 7 uses gaps/shape_and_iteration_tables/barrido_FK0_TE.dat.
  Fig. 8 uses the stored equivalent-area factors and sensitivity
  bounds in gaps/shape_and_iteration_tables/ahorro_consistente_TE.csv.
  Fig. S2 uses the stored IEM and FEM differences from
  validation/fem_pares_TE.dat.

Figs. 9 and 11
  Command  python3 figures/plot_summary_figures.py OUTPUT
  Fig. 9 uses the reference crossing, sensitivity limits, and paired
  refinements in gaps/saving_methods.csv, gaps/smooth_curves.csv,
  gaps/paired_factors.csv and gaps/paired_curves.csv.
  Fig. 11 uses fields/recess_by_mode.csv and fields/recess_refinement.csv.

Fig. 10 and Figs. S3–S5
  Command  python3 figures/plot_field_figures.py OUTPUT
  These figures use the fields/manuscript_maps.npz intensity arrays,
  the associated figure index and metadata, stored widths in
  fields/map_bandwidths_and_controls.csv, and stored overlaps in
  fields/intensity_correlations.json.
  Field maps are exported at 500 dpi in the separate PDFs. Labels,
  contours and plotted profiles retain their vector representation.

Fig. 12 and the geometry-block resampling
  Command  python3 bands/plot_coupling.py OUTPUT/bandwidth_and_channel_attenuation.pdf
  This script draws the fitted lines from bands/triangle_modes.csv
  and bands/diagnostics/acoplamiento_wt_TE.dat using the stated fits.
  The optional flag --check-bootstrap repeats the resampling check
  (10,000 draws, full geometries as blocks, seed 0). Each selected
  geometry carries all seven of its modes.

Integrity
MANIFEST.sha256 lists the SHA-256 checksum of each included file.
Run shasum -a 256 -c MANIFEST.sha256 from the dataset root to verify it.
