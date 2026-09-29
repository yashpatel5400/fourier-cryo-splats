"""Generate manuscript results and a human-readable report solely from saved runs."""

from pathlib import Path
import json
import re
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
datasets = ["10028", "10049", "10076"]
metrics = {
    d: json.loads((ROOT / "results/final" / d / "metrics.json").read_text())
    for d in datasets
}
baselines = {
    d: json.loads((ROOT / "results/final" / d / "cryodrgn-baseline.json").read_text())
    for d in datasets
}
provenance = {
    d: json.loads((ROOT / "provenance" / (d + "-manifest.json")).read_text())
    for d in datasets
}
noise = {
    d: json.loads((ROOT / "results/noise" / d / "metrics.json").read_text())
    for d in datasets
}
selection = json.loads((ROOT / "results/selection.json").read_text())
synthetic = json.loads((ROOT / "results/validation/synthetic.json").read_text())
adaptive = json.loads((ROOT / "results/adaptive-demo/10049/metrics.json").read_text())
truncation = json.loads((ROOT / "results/validation/truncation.json").read_text())


def res(r, tex=False):
    if r["angstrom"] is None:
        return "--"
    prefix = (
        ("\\leq " if tex else "≤") if r["status"] == "censored_at_sampled_limit" else ""
    )
    value = prefix + f"{r['angstrom']:.2f}"
    return "$" + value + "$" if tex else value


lines = [
    r"""\subsection{Executed study}
We downloaded real deposited particle bytes for 8,192 particles from each of
three distinct accessions, totaling 24,576 experimental images. Selection used
randomly permuted contiguous blocks of 64 in the source order, applying the
published cryoDRGN particle filter when present. This is a reproducible block
sample, not an independent uniform sample or the entire deposited dataset.
Original selected images were Fourier-cropped to $64\times64$ while preserving
the physical field of view. Source byte ranges, SHA-256 hashes, particle indices,
and the metadata repository commit are included in the reproducibility package.
No detector-movie motion correction or new particle picking was performed.
One EMPIAR-10028 source file was excluded because its header declared 113
particles but the served file contained only 4 KB. The 32 affected selected
particles were replaced by continuing the same seeded selection order; the
exclusion and identity-preserving recovery are documented in the repository.

\begin{table}[t]
\caption{Experimental inputs. All runs use 8,192 particles per accession. Pixel
sizes are in \AA/pixel; ``raw'' denotes deposited extracted images.}
\label{tab:data}
\centering\small
\begin{tabular}{lrrr}\toprule
EMPIAR & Raw box & Raw pixel & Working pixel\\\midrule"""
]
for d in datasets:
    m = provenance[d]
    lines.append(
        f"{d} & {m['raw_box']} & {m['raw_pixel_size_A']:.2f} & {m['pixel_size_A']:.4f} \\"
    )
lines += [
    r"""\bottomrule\end{tabular}\end{table}
For each accession, a fixed permutation (seed 42 for observations; 43 for
partitioning) reserves 409 images for validation and 409 for test. The remaining
7,374 images are split into two disjoint sets of 3,687. We use all 1,410 nonredundant
Fourier samples inside radius 30, excluding DC. Images are divided by their
background standard deviation and given the documented cryoDRGN sign convention.
The standard linear image-window taper from 0.85 to 0.99 of the half-box is used
for the final experiment. Windowing introduces Fourier noise correlations; our
diagonal least-squares objective is therefore an approximation, shared by the
matched methods. The final map is the arithmetic mean of the two fitted maps.

The stationary Gaussian dictionary has centers on a $65^3$ lattice, conjugate
tied across the origin, width 0.5 Fourier bins, and radial truncation at 2 bins.
The implementation allocates 274,625 real scalar coefficient slots, including
unsupported locations driven to zero. It does not learn these centers or widths
in the main three-dataset experiment. Preconditioned conjugate gradients solve
the real and imaginary systems separately with relative tolerance $10^{-5}$ and
at most 80 iterations. The ridge coefficient is $\lambda=\alpha$
times the median positive diagonal of the corresponding data normal matrix.

The matched baseline replaces the Gaussian basis by trilinear Fourier voxels
while preserving data, frequency samples, objective, partitions and solver.
A second baseline invokes \texttt{cryoDRGN 4.3.1 backproject\_voxel}, with
CTF multiplication, regularization weight 1, the same particle halves and the
software's standard image window. That external baseline uses its native
backprojection rule and circular frequency support; FSC comparisons restrict all
maps to the same radius-30 sphere. It is a classical fixed-pose command, not a
neural or ab initio cryoDRGN experiment.

\subsection{Regularization selected by prediction}
We selected $\alpha$ from $\{0.1,1,10\}$ using RAG validation-image error,
independently for each representation. Both selected $\alpha=1$, which was then
frozen for the other accessions. Neither the reported FSC curves nor test error
selected this value. Table~\ref{tab:reg} shows why a high FSC alone is insufficient:
the strongest regularizer can preserve high correlation while underpredicting
the observed signal.
\begin{table}[t]
\caption{RAG validation NMSE used for selecting regularization. Lower is better.}
\label{tab:reg}\centering\small
\begin{tabular}{rrr}\toprule
$\alpha$ & Gaussian & Matched voxel\\\midrule"""
]
for a in ["0.1", "1", "10"]:
    lines.append(
        f"{a} & {selection['candidates']['gaussian'][a]:.6f} & {selection['candidates']['voxel'][a]:.6f} \\"
    )
lines += [
    r"""\bottomrule\end{tabular}\end{table}

\subsection{Reconstruction agreement and FSC}
Figure~\ref{fig:fsc} reports unmasked conditional half-map FSC. No spatial map
mask, sharpening, or test-derived alignment is used. Table~\ref{tab:fsc} distinguishes
a measured threshold crossing from right-censoring at the last sampled frequency.
A censored value is a bound at this experiment's sampling limit, not a measurement
of a finer resolution. Figure~\ref{fig:maps} shows representative map sections.

\begin{figure*}[t]\centering
\includegraphics[width=\textwidth]{../figures/final-fsc.png}
\caption{Conditional half-map FSC on three experimental particle selections.
Poses are supplied from public consensus refinements, so these are not fully
independent gold-standard reconstructions. The dashed line is 0.143.}
\label{fig:fsc}\end{figure*}

\begin{table}[t]
\caption{Conditional half-map FSC at 0.143, in \AA. A $\leq$ sign means the
curve never crosses within the sampled band. BP is cryoDRGN's classical
backprojection. No result in this table is an ab initio resolution estimate.}
\label{tab:fsc}\centering\small
\begin{tabular}{lrrr}\toprule
EMPIAR & Gaussian & Voxel & BP\\\midrule"""
]
for d in datasets:
    j = metrics[d]
    lines.append(
        f"{d} & {res(j['methods']['gaussian']['half_fsc_resolution'], True)} & {res(j['methods']['voxel']['half_fsc_resolution'], True)} & {res(baselines[d]['half_fsc_resolution'], True)} \\"
    )
lines += [
    r"""\bottomrule\end{tabular}\end{table}

\begin{table}[t]
\caption{Mean shell-wise cross-method FSC over shells 1--30. These are agreement
statistics between methods sharing images and poses, not independent accuracy
or resolution measurements.}
\label{tab:agreement}\centering\small
\begin{tabular}{lrr}\toprule
EMPIAR & Gaussian--voxel & Gaussian--BP\\\midrule"""
]
for d in datasets:
    lines.append(
        f"{d} & {metrics[d]['mean_method_agreement_fsc']:.4f} & {baselines[d]['cross_method_mean_fsc']:.4f} \\"
    )
lines += [
    r"""\bottomrule\end{tabular}\end{table}

\begin{figure*}[t]\centering
\includegraphics[width=\textwidth]{../figures/final-slices.png}
\caption{Central sections through the reconstructed mean maps. Rows are
EMPIAR-10028, 10049, and 10076; columns show three orthogonal sections of Gaussian
and matched voxel maps. For visualization only, each volume is smoothed by one
voxel and displayed on its own intensity scale. FSC uses the unsmoothed maps.}
\label{fig:maps}\end{figure*}

\subsection{Held-out prediction and cost}
Table~\ref{tab:prediction} reports normalized mean squared error on the 409 test
images excluded from fitting and regularizer selection. Values close to one are
expected for noisy particle images, since the denominator includes measurement
noise. A small improvement should not be interpreted as a large structural
advance. Supplemental particle-bootstrap intervals are conditional on the
supplied poses and do not account for correlations within micrographs.

\begin{table}[t]
\caption{Held-out test NMSE and elapsed reconstruction time. G/V denote Gaussian
and matched voxel. Timings include sparse operator construction, solves, map
export, and held-out prediction; they exclude downloads and initial data loading.}
\label{tab:prediction}\centering\small
\begin{tabular}{lrrrr}\toprule
EMPIAR & G NMSE & V NMSE & G (s) & V (s)\\\midrule"""
]
for d in datasets:
    g = metrics[d]["methods"]["gaussian"]
    v = metrics[d]["methods"]["voxel"]
    lines.append(
        f"{d} & {g['test_nmse']:.5f} & {v['test_nmse']:.5f} & {g['seconds']:.1f} & {v['seconds']:.1f} \\"
    )
lines += [
    r"""\bottomrule\end{tabular}\end{table}
The main runs use CPU sparse linear algebra on an Apple M4 Pro with 24 GiB
unified memory. The reference Gaussian implementation is slower than the
trilinear baseline, with much of its cost in sparse kernel construction.
These runs do not establish a speed or memory advantage. Concurrent downloads
and operating-system activity also limit the precision of timing comparisons.

\subsection{Controls and a nonlinear feasibility run}
An out-of-dictionary synthetic phantom is generated from displaced real-space
Gaussians. The low-shell mean FSC of the recovered Fourier Gaussian model is
""",
    f"{synthetic['gaussian']['low_shell_mean_fsc']:.4f}, versus {synthetic['voxel']['low_shell_mean_fsc']:.4f} for the voxel baseline.",
    r"""
This is an implementation check, not an experimental-data accuracy result.
Analytic plane restriction, Hermitian symmetry, sparse adjoints, CTF parity,
FSC behavior, inverse-transform reality and nonlinear geometry gradients also
pass numerical checks. Increasing the Gaussian evaluation cutoff from 2 to 3
bins at fixed fitted coefficients produces the relative complex prediction
errors recorded with the released validation outputs; this checks rendering
truncation, not the effect of refitting with a different basis.

A phase-randomized control is fitted independently on each accession at
$32\times32$, with 256 Fourier samples per image. This preserves Fourier
magnitudes while destroying coherent phase information. We report its entire
FSC curves, rather than a single favorable shell. None of these controls is a
substitute for independent half-set pose refinement.

The general differentiable representation is separately tested on low-frequency
RAG observations using a $16\times16$ grid, """,
    str(adaptive["gaussian_pairs"]),
    r""" conjugate pairs, and
joint optimization of complex coefficients, centers and full anisotropic
precisions. This MPS run uses 7,680 training images, 512 validation images, 1,500
Adam steps, bounded center displacements, and bounded precision parameters.
The initial and best validation NMSE are """,
    f"{adaptive['initial_validation_nmse']:.4f} and {adaptive['best_validation_nmse']:.4f}",
    r""", respectively,
with """,
    f"{adaptive['seconds']:.1f}",
    r""" seconds of optimization. Its partition and band differ
from the main benchmark. This demonstrates nonlinear fitting feasibility only;
it does not establish improved resolution, sparse adaptive efficiency, or
superiority to the fixed dictionary.""",
]
# Each table row must use the LaTeX double backslash, not a single line escape.
text = "\n".join(lines)
figure_pattern = r"\\begin\{figure\*\}.*?\\end\{figure\*\}"
figure_blocks = re.findall(figure_pattern, text, flags=re.S)
text = re.sub(figure_pattern, "", text, flags=re.S)
# Declare double-column floats early so they appear alongside results, before references.
text = "\n".join(figure_blocks) + "\n" + text
max_truncation = max(v["relative_complex_L2_error"] for v in truncation.values())
noise_values = []
for dataset in datasets:
    curve = np.loadtxt(
        ROOT / "results/noise" / dataset / "gaussian-half-fsc.csv",
        delimiter=",",
        skiprows=1,
    )
    noise_values.append(float(np.mean(np.abs(curve[3:, 2]))))
text = text.replace(
    "This checks rendering truncation", "This checks rendering truncation"
)
text = text.replace(
    "recorded with the released validation outputs",
    f"recorded with the released validation outputs (maximum {max_truncation:.5f} across datasets)",
)
text = text.replace(
    "None of these controls is a",
    "The mean absolute Gaussian-control FSC over shells 4--14 is "
    + ", ".join(f"{v:.3f}" for v in noise_values)
    + " for the three accessions, respectively. None of these controls is a",
)
text = "\n".join(
    line + "\\" if line.endswith(" \\") and not line.endswith(" \\\\") else line
    for line in text.split("\n")
)
text = "\n".join(line.rstrip() for line in text.splitlines())
(ROOT / "paper/results.tex").write_text(text + "\n")
summary = []
for d in datasets:
    j = metrics[d]
    summary.append(
        {
            "dataset": d,
            "particles": j["particle_count"],
            "gaussian_half_fsc": j["methods"]["gaussian"]["half_fsc_resolution"],
            "voxel_half_fsc": j["methods"]["voxel"]["half_fsc_resolution"],
            "cryodrgn_backprojection_half_fsc": baselines[d]["half_fsc_resolution"],
            "gaussian_vs_voxel_mean_fsc": j["mean_method_agreement_fsc"],
            "gaussian_vs_backprojection_mean_fsc": baselines[d][
                "cross_method_mean_fsc"
            ],
            "gaussian_test_nmse": j["methods"]["gaussian"]["test_nmse"],
            "voxel_test_nmse": j["methods"]["voxel"]["test_nmse"],
            "gaussian_seconds": j["methods"]["gaussian"]["seconds"],
            "voxel_seconds": j["methods"]["voxel"]["seconds"],
        }
    )
(ROOT / "results/summary.json").write_text(json.dumps(summary, indent=2))
md = [
    "# Measured results",
    "",
    "All half-map FSC values are conditional on public consensus poses. A ≤ sign denotes a sampling-limit bound, not a measured finer resolution.",
    "",
    "| EMPIAR | Gaussian FSC (Å) | Matched voxel (Å) | cryoDRGN BP (Å) | Gaussian vs BP mean FSC |",
    "|---|---:|---:|---:|---:|",
]
for d in datasets:
    md.append(
        f"| {d} | {res(metrics[d]['methods']['gaussian']['half_fsc_resolution'])} | {res(metrics[d]['methods']['voxel']['half_fsc_resolution'])} | {res(baselines[d]['half_fsc_resolution'])} | {baselines[d]['cross_method_mean_fsc']:.4f} |"
    )
md += [
    "",
    "24,576 real particles total; 8,192 per accession; 64×64 images. The model uses fixed Gaussian geometry in the main experiment. It is slower than the matched voxel solver. No ab initio, full-stack, near-atomic, heterogeneity, or fully gold-standard claim is made.",
    "",
    "The public paper and JSON/CSV files contain the full protocol, predictive error, timing, ablations, and limitations.",
]
(ROOT / "RESULTS.md").write_text("\n".join(md) + "\n")
print(json.dumps(summary, indent=2))
