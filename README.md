# Endophyte-dependent transcriptional responses to a pathogen in Norway spruce

This repository accompanies the manuscript *Endophyte-dependent transcriptional responses to a pathogen in Norway spruce*. It contains the derived data tables, analysis scripts, sequence-domain inputs and models, and provenance records used for the reported two-factor transcriptomic reanalysis of BioProject PRJNA734471.

## Scope

The archived materials support inspection and reproduction of the reported derived analyses. The original RNA-seq reads are not redistributed here. They remain available through NCBI SRA BioProject PRJNA734471. Large reference resources and their retrieval records are described in the included provenance files.

## Repository layout

- `data/`: electronic Supplementary Table S5 and derived candidate-level data.
- `analysis/`: analysis results, scripts, enrichment records, domain-verification materials, and receipts.
- `frozen_base/`: frozen input, provenance, numerical, and reproducibility records.
- `reproducibility_inputs/`: archived model objects and inputs for the derived analyses.
- `figures/`: figure-related outputs included in the original reproducibility package.
- `work/`: supporting scripts used to produce or check archived derived outputs.

Electronic Supplementary Table S5 is provided in editable tab-separated format: [candidate-level results](data/S5_candidates39_complete.tsv) and [source GO relations](data/S5_source_GO_relations.tsv). See the [S5 readme](data/README_ELECTRONIC_TABLE_S5_EN.txt) for field descriptions.

## Reproduction

The archived scripts retain the original analytical record and may contain historical absolute paths. Set local paths before executing them. Do not overwrite the supplied frozen results. Core workflows use R (including DESeq2 and jsonlite); Figure 3 additionally uses ggplot2, data.table, and grid. The domain-verification workflow uses Python with pyhmmer. See `analysis/README_REPRODUCE.txt` and `README_REPRODUCE_EN.txt` for the detailed workflow and input provenance.

## Data and provenance

The repository redistributes only files assembled for this study's reproducibility package. It does not grant rights to third-party sequence reads, database models, annotations, or external reference resources. Consult the source-accession and provenance records before retrieving or reusing those materials.

## Citation

Please cite the associated manuscript and this archived repository version: Ma, Z., and Liu, Z. (2026). *Endophyte-dependent transcriptional responses to a pathogen in Norway spruce: code and derived data*. Version 1.0.0. https://github.com/mzw423502/dse-spruce-transcriptome-reanalysis/tree/v1.0.0 . A DOI has not been assigned; use `CITATION.cff` for citation metadata.
