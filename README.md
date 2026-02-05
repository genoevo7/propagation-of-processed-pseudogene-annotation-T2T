# Propagation of Processed Pseudogene Annotation in the T2T Era

Code repository accompanying the paper:

Boundary-associated propagation of a processed pseudogene dissects pre-existing limitations of genome annotation in the T2T era  
Published in Mobile DNA (2026)

---

## Overview

This repository provides all Python scripts used to reproduce the computational analyses and figures presented in the Mobile DNA (2026) paper:

Boundary-associated propagation of a processed pseudogene dissects pre-existing limitations of genome annotation in the T2T era.

The study investigates how a processed pseudogene locus propagates in association with chromosomal boundaries in the telomere-to-telomere (T2T) genome era, and how such propagation exposes pre-existing limitations in genome annotation frameworks.

The repository focuses exclusively on reproducible, block-based, and projection-aware analyses. All external inference steps (e.g., selection tests) are explicitly separated from in-repository visualization and quantification.

---

## Repository Contents

blocks_locus_dispersion.py  
blocks_chr_painting.py  
gapmasked_alignment_fel_meme_plot.py  
requirements.txt  
README.md  

---

## 1. Block-Based Locus Dispersion Analysis

Script: blocks_locus_dispersion.py

This script quantifies how a defined focal locus in the reference genome is distributed across query genomes using UCSC chain/net block projections.

The analysis is performed in a genome-wide, projection-based framework and does not assume true insertion events. All results are interpreted strictly as projected block mappings.

### Inputs

UCSC Table Browser chain/net block exports (.txt), one file per query species.

Expected columns:
chrom, chromStart, chromEnd, strand, qName, qSize, qStart, qEnd

Both header and no-header formats are supported.

### Target locus definition

The focal locus is defined directly in the script, for example:

TARGET_TNAME  = chr7  
TARGET_TSTART = 55953761  
TARGET_TEND   = 55959097  

### Outputs

all_species_locus_by_qname.tsv  
all_species_locus_species.tsv  

### Reported metrics

- dispersion_qnames: number of distinct query chromosomes or contigs
- total_hit_bp: cumulative projected block length
- dupIndex: total_hit_bp / window_len_bp
- mix_fraction: overlap between plus/minus strand blocks
- conflict_fraction: strand conflict normalized by union length

These metrics are used to quantify locus dispersion and structural amplification without asserting functional or insertional causality.

---

## 2. Chromosome Painting Visualization

Script: blocks_chr_painting.py

This script generates chromosome painting plots showing how blocks derived from a focal human locus are distributed across chromosomes or contigs in each query species.

The script robustly handles assemblies containing haplotypes or patch sequences. When chromosome labels (e.g., chr7) are duplicated across haplotypes (e.g., hap1 / hap2), the script automatically switches to an internal chromosome-label keyspace to prevent name collisions and ensure successful plotting.

### Inputs

- UCSC chain/net block files (same as above)
- Chromosome size tables downloaded from NCBI Datasets

Chromosome size tables contain chromosome (or contig) names and lengths.

### Gene coordinate annotation

Gene coordinates are manually curated from the NCBI Genome Browser for each species and assembly and explicitly defined inside the script.

Example:

GENE_LOCS = {
  "Chimpanzee chromosome painting": {
    "SEPTIN14": ("chr6", 61565981, 61648745),
    "CIC":      ("chr20", 45184628, 45211932),
  }
}

These coordinates are entered manually to ensure assembly-consistent annotation and are not inferred automatically.

### Outputs

High-resolution chromosome painting figures (PNG and PDF), suitable for direct inclusion in publication figures.

---

## 3. Alignment-Aware Visualization of Codon-Level Selection Signals

Script: gapmasked_alignment_fel_meme_plot.py

This script visualizes codon-aligned CDS regions across multiple species, integrating alignment gaps and externally inferred selection signals.

The visualization maps selection signals onto an alignment-aware CDS framework and does not perform selection inference itself.

### External workflow (performed outside this repository)

1. Ortholog CDS retrieval  
   Representative species are selected using NCBI Ortholog resources. CDS sequences are downloaded, and stop codons are removed.

2. Codon-based alignment  
   Codon-based MUSCLE alignment is performed using MEGA12. The resulting alignment is exported as a FASTA file.

3. Selection inference  
   The codon alignment is uploaded to the Datamonkey web server, and selection analyses are performed using HyPhy:
   - FEL (Fixed Effects Likelihood)
   - MEME (Mixed Effects Model of Evolution)

The resulting Excel files are used as inputs for visualization in this repository.

### Inputs

- Codon-based CDS alignment FASTA file (stop codons removed)
- FEL result Excel file
- MEME result Excel file

### Outputs

Gap-masked, alignment-aware figures summarizing codon-level selection signals across species, exported as PNG and PDF files.

---

## Dependencies

The following Python packages are required to run the scripts in this repository:

numpy>=1.21  
pandas>=1.4  
matplotlib>=3.5  
biopython>=1.79  
openpyxl>=3.0  

Install dependencies with:

pip install -r requirements.txt

---

## Reproducibility and Interpretation Notes

- All analyses are based on UCSC chain/net projections.
- No claim is made that projected blocks represent true insertion events.
- Selection inference is performed externally; this repository is limited to visualization and structural interpretation.
- All coordinates, parameters, and assumptions are explicitly defined in the scripts.

This design ensures transparency and avoids over-interpretation of projection-based signals.

---

## License

This repository is provided for academic and research use.  
If you use this code or figures in published work, please cite the corresponding Mobile DNA (2026) article.
