# Study results

Keep compact summaries, statistics, manifests, reports and curated plot inputs
in Git. The bulky CARMA per-run CSV tables are generated, ignored build outputs.
A report-only stage needs those local tables: rerun its producing stages first
on a fresh checkout. Ignoring them does not disable writing or resuming studies.

Build `tests/kalman_carma/carma3_bench` with
`make -C tests/kalman_carma carma3_bench` before the CARMA studies. The shared-mode
studies also require `make -C tests/kalman_carma_shared_modes build`. Use the
seed counts and options recorded in each study's manifest and research document
when reproducing published results; driver defaults can differ.

| Directory | Driver under `tools/` | Stages producing ignored tables |
| --- | --- | --- |
| `carma3_benchmark` | `carma3_benchmark.py` | `run` |
| `carma3_coupled` | `carma3_coupled_benchmark.py` | `screen`, `run`, `ablate` |
| `carma_excitation` | `carma_excitation_study.py` | `tune`, `run` |
| `tfg_carma3` | `tfg_carma3_benchmark.py` | `run`, `followup` |
| `kalman_carma_shared_modes` | `carma_shared_modes_validation.py` | `screen`, `run`, `robust` |
| `kalman_carma_shared_modes_minimal` | `carma_shared_modes_minimal_study.py` | `run`, `ablate`, `robust` |

For example, `python3 tools/carma3_benchmark.py run` recreates `raw_runs.csv`;
`python3 tools/carma3_benchmark.py traces` recreates the local figure traces;
`python3 tools/carma3_benchmark.py report` builds the report and SVG figures.
Other stages such as dataset replay, spectra and timing have separate inputs and
prerequisites documented in the drivers. Manifests describe generated runs and
may list local outputs that are intentionally absent from Git.

The OU validation/robustness JSON bundles and transition CSVs remain versioned:
report rebuilding and provenance checks consume them. File size alone does not
make these artifacts redundant.

All Matplotlib exporters in `plots/` and `tools/` use
`plots/vector_export.py`. Each dense line, 3D trajectory or filled-band polygon
is reduced to at most 2,000 vertices for export, retaining ordered bucket
coordinate extrema, endpoints and nonfinite breaks. SVG, PGF and PDF receive the
same reduced geometry, so downstream LaTeX PDFs do not embed full time series.
Dense overlaid curves share the axes budget, with at least 256 points per
curve. Short curves are unchanged. Sampling is visual only: input CSVs, FFTs, scores,
axis limits and the in-memory source artists retain their original data.
The existing rate-based input sampling can still reduce load time; the export
budget additionally bounds charts whose duration grows on later study reruns.

Removing generated files from the current tree does not shrink Git history.
