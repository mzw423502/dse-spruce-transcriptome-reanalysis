SUPPORTING DATA AND CODE: REPRODUCIBILITY GUIDE

Scope
-----
This archive contains the frozen inputs, statistical exports, functional-source
ledgers, targeted Pfam-domain materials, GO biological-process inputs, scripts,
and final tabular outputs used for the transcriptomic reanalysis of Norway
spruce responses to dark septate endophyte (DSE) and pathogen treatments. It is
designed to permit transparent reconstruction and inspection of the reported
derived tables and figures from frozen products. It does not redistribute raw
FASTQ files, a preconfigured software environment, or manuscript and editorial
records. Raw sequencing reads are publicly available from NCBI SRA BioProject
PRJNA734471 and should be retrieved from the public repository when raw-read
processing is required.

Archive layout
--------------
data/
  Electronic Table S5 source files. The two TSV files are retained byte-for-byte
  from the frozen delivery.
frozen_base/
  Minimal frozen numerical and formal-statistics inputs required by the replay
  scripts, plus the functional-source ledgers required for candidate provenance.
reproducibility_inputs/
  The already fitted primary DESeq2 model object used for descriptive normalized
  count extraction. This archive does not refit the RNA-seq model.
work/
  Entry-point scripts for classification/expression replay, source-evidence
  extraction, candidate-table assembly, and Figure 1/2 rendering.
analysis/
  Frozen targeted Pfam-domain inputs/models/results and GO-BP input, scripts,
  result tables, validation tables, and the portable Figure 3 renderer.

Software requirements
---------------------
R (with DESeq2 and jsonlite) is required for the core replay. Base R is
sufficient for the Figure 1/2 renderer. GO-BP analysis and Figure 3 rendering
also require R packages ggplot2, data.table, and grid. Python 3 is required for
the source audit and validation scripts; numpy is required by the GO validator.
The Pfam rerun scripts require pyhmmer and the included HMM files. These
dependencies are not bundled in the archive.

Core replay
-----------
Create new, empty output directories outside the archive. From the archive
root, run:

  mkdir replay_analysis replay_review replay_figures
  Rscript --vanilla work/reproduce_revision.R frozen_base reproducibility_inputs/primary_DESeq2_model.rds replay_analysis
  python work/reproduce_source_audit.py frozen_base replay_review
  python work/reproduce_candidate_table.py replay_analysis replay_review/opposite39_source_function_audit.tsv replay_review/S5_candidates39_replayed.tsv
  Rscript --vanilla work/reproduce_figures.R replay_analysis replay_figures

The first command reconstructs derived classifications and descriptive
expression exports from frozen result tables. It does not quantify reads or
refit the primary model. The second and third commands recreate the
source-evidence join and the candidate table. The figure script consumes the
reconstructed tables and does not perform an inferential reanalysis.

GO-BP and Figure 3
------------------
Copy analysis/enrichment to a writable working location before running its
validator. Run the GO script with a new output directory; for example:

  Rscript --vanilla analysis/enrichment/scripts/run_go_statistics.R --input-dir analysis/enrichment/input --output-dir replay_go --n-perm 10000 --seed 20261110
  cp -r analysis/enrichment replay_enrichment
  python replay_enrichment/scripts/validate_go_results.py replay_enrichment
  Rscript --vanilla analysis/figure_revision/plot_figure3_enrichment_candidates_portable.R analysis/enrichment replay_figure3

The validator accepts the GO project root and writes its validation products
under that project root. The Figure 3 command reads the frozen result and
validation tables in analysis/enrichment/run. The GO rerun creates a separate
result directory and does not alter the supplied frozen outputs.

Targeted Pfam-domain rerun
--------------------------
After installing pyhmmer in the active Python environment, run:

  python analysis/rerun_from_frozen_inputs.py --root analysis --output replay_domain
  python analysis/rerun_diagnostic_from_frozen_inputs.py --root analysis --output replay_domain_diagnostic

The first command applies the included Pfam HMMs at their predefined gathering
thresholds. The diagnostic command reports the separately designated outside-GA
records. Neither command changes the frozen tables.

Path and execution notes
------------------------
The listed entry-point scripts use explicit archive-relative arguments. Some
legacy source scripts from the original computational workspace used absolute
paths; they are not execution entry points in this public archive. If a retained
script is adapted locally and contains an absolute path, replace that path with
the local archive root without changing frozen data values. This package is a
reproducibility resource rather than a fully preconfigured one-command software
environment.

Data boundaries
---------------
The analysis concerns a single public cohort. The supplied sensitivity analyses
evaluate internal consistency across analytical routes; they do not constitute
independent biological replication. Functional labels and targeted domain
results should be interpreted according to the definitions in the manuscript
and supplementary information.
