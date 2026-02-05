# Locus Dispersion, Chromosome Painting, and Selection Analysis Pipelines

This repository contains Python scripts used for block-based locus dispersion quantification, chromosome painting visualization, and codon-level selection result visualization across multiple primate species. All analyses are designed to be transparent, reproducible, and explicitly separated from external (non-Python) tools.

---

## Overview

This repository provides three independent but complementary analysis components:

1. Block-based locus dispersion analysis  
   Quantifies how a defined human genomic locus is dispersed across query genomes using UCSC chain/net block data.

2. Chromosome painting visualization  
   Visualizes how a focal locus maps onto multiple chromosomes or contigs in query species, including cases with haplotypes or patch assemblies.

3. Codon-level selection visualization (FEL / MEME)  
   Visualizes purifying and episodic diversifying selection sites inferred externally using Datamonkey (HyPhy), mapped onto codon or CDS coordinates.

External steps such as codon-based alignment and selection inference are explicitly separated from Python-based visualization.

---

## Repository Contents

blocks_locus_dispersion.py  
blocks_chr_painting.py  
selection_sites_FEL_MEME.py  
Figure_SEPTIN14_gapmasked_FEL_MEME.py  
requirements.txt  
README.md  

---

## 1. Block-Based Locus Dispersion Analysis

Script: blocks_locus_dispersion.py

This script quantifies how a defined target locus in the reference genome (t) is distributed across query genomes (q) using UCSC Table Browser chain/net block exports.

Metrics include:
- Number of distinct query chromosomes or contigs (dispersion_qnames)
- Total aligned length in query genomes
- Strand conflicts between plus and minus blocks
- Duplication / dispersion index

Input files are UCSC Table Browser exports (.txt), placed in the same directory. Header and no-header formats are both supported.

Expected columns:
chrom, chromStart, chromEnd, strand, qName, qSize, qStart, qEnd

Example input files:
chimpanzee.txt  
pygmy chimpanzee.txt  
western lowland gorilla.txt  
Bornean orangutan.txt  
Sumatran orangutan.txt  
siamang.txt  
slow loris.txt  
Ring-tailed lemur.txt  
white-tufted-ear marmoset.txt  

The target locus is defined directly in the script:

TARGET_TNAME  = chr7  
TARGET_TSTART = 55953761  
TARGET_TEND   = 55959097  

Outputs:
all_species_locus_by_qname.tsv  
all_species_locus_species.tsv  

Key metrics:
dupIndex = total_hit_bp / window_len_bp  
mix_fraction = overlap(+,-)_bp / hit_bp  
conflict_fraction = overlap(+,-)_bp / union(+,-)_bp  

---

## 2. Chromosome Painting Visualization

Script: blocks_chr_painting.py

This script produces chromosome painting plots showing how blocks derived from a focal human locus are distributed across chromosomes or contigs in each query species.

The script automatically handles cases where chromosome labels (e.g., chr7) are duplicated across haplotypes or patch assemblies (e.g., hap1 / hap2) by switching internally to a chromosome-label keyspace, ensuring that plots are generated robustly without label collisions.

Inputs:
- UCSC block files (same as above)
- Chromosome size tables downloaded from NCBI Datasets

Gene coordinates are manually curated from the NCBI Genome Browser for each species and assembly and defined explicitly inside the script.

Example:

GENE_LOCS = {
  "Chimpanzee chromosome painting": {
    "SEPTIN14": ("chr6", 61565981, 61648745),
    "CIC":      ("chr20", 45184628, 45211932),
  }
}

These coordinates are manually inspected and entered based on the corresponding NCBI reference assemblies.

Outputs:
High-resolution chromosome painting figures (PNG and PDF), suitable for publication.

---

## 3. Codon-Level Selection Visualization (FEL / MEME)

Selection inference itself is not performed in this repository. Instead, this repository visualizes results generated using established external tools.

External workflow:

1. Ortholog CDS retrieval  
   Representative species are selected from NCBI Ortholog. CDS sequences are downloaded, and stop codons are removed.

2. Codon-based alignment  
   Codon-based MUSCLE alignment is performed using MEGA12. The resulting alignment is exported as a FASTA file.

3. Selection inference  
   The codon alignment is uploaded to Datamonkey, and selection analyses are performed using HyPhy:
   - FEL (Fixed Effects Likelihood)
   - MEME (Mixed Effects Model of Evolution)

The resulting Excel (.xlsx) files are used as inputs for Python-based visualization.

---

## 3a. FEL / MEME Codon-Level Plot

Script: selection_sites_FEL_MEME.py

Inputs:
- FEL result Excel file
- MEME result Excel file

Expected columns:

FEL:
codon, p-value, class

MEME:
Codon, p-value, Class

Outputs:
selection_sites_FEL_MEME.png  
selection_sites_FEL_MEME.pdf  

This script produces Manhattan-style codon-level plots with capped −log10(p) values and optional highlighted positions.

---

## 3b. Alignment-Aware Gap-Masked CDS Visualization

Script: Figure_SEPTIN14_gapmasked_FEL_MEME.py

This script visualizes codon-aligned CDS blocks with gap-masked regions, FEL purifying sites, and MEME episodic diversifying sites across multiple species.

Inputs:
- Codon-based alignment FASTA file (stop codons removed)
- FEL result Excel file
- MEME result Excel file

Outputs:
Figure_SEPTIN14_gapmasked_FEL_MEME.png  
Figure_SEPTIN14_gapmasked_FEL_MEME.pdf  

---

## Dependencies

Python dependencies required only for running the scripts in this repository:

numpy>=1.21  
pandas>=1.4  
matplotlib>=3.5  
biopython>=1.79  
openpyxl>=3.0  

Install with:
pip install -r requirements.txt

---

## Notes on Reproducibility

All Python scripts are deterministic given the same inputs. External tools (MEGA12 and Datamonkey/HyPhy) are explicitly documented but not bundled. No hidden preprocessing steps are performed. All coordinates, thresholds, and assumptions are visible in the code.

---

## License

This code is provided for academic and research use. Please cite appropriately if reused in published work.
