# DSE-dependent transcriptomic responses to a pathogen in Norway spruce

This repository accompanies the manuscript *Endophyte-dependent transcriptional responses to a pathogen in Norway spruce*. It contains the derived data tables, analysis scripts, sequence-domain inputs and models, and provenance records used for the reported two-factor transcriptomic reanalysis of BioProject PRJNA734471.

## Scope

The archived materials support inspection and reproduction of the reported derived analyses. The original RNA-seq reads are not redistributed here. They remain available through NCBI SRA BioProject PRJNA734471. Large reference resources and their retrieval records are described in the included provenance files.

## Repository layout

- `data/`: electronic Supplementary Table S5 and derived candidate-level data.
- `analysis/`: analysis results, scripts, enrichment records, domain-verification materials, and receipts.
- `frozen_base/`: frozen input, provenance, numerical, and reproducibility records.
- `work/`: supporting scripts used to produce or check archived derived outputs.

## Reproduction

The archived scripts retain the original analytical record and may contain historical absolute paths. Set local paths before executing them. Do not overwrite the supplied frozen results. Core workflows use R (including DESeq2 and jsonlite); Figure 3 additionally uses ggplot2, data.table, and grid. The domain-verification workflow uses Python with pyhmmer. See `analysis/README_REPRODUCE.txt` and `README_REPRODUCE_EN.txt` for the detailed workflow and input provenance.

## Data and provenance

The repository redistributes only files assembled for this study's reproducibility package. It does not grant rights to third-party sequence reads, database models, annotations, or external reference resources. Consult the source-accession and provenance records before retrieving or reusing those materials.

## Citation

Please cite the associated manuscript. A DOI is not assigned at the time of repository preparation; use `CITATION.cff` for the current citation metadata.
