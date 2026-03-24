# Hydrophobicity Space

A toolkit for defining and visualizing amino acid hydrophobicity across 23 published hydrophobicity scales, scraped from [ExPASy ProtScale](https://web.expasy.org/protscale/).

## Overview

Different hydrophobicity scales assign different values to amino acids depending on the experimental method used (partitioning, chromatography, structural analysis, etc.). This project aggregates 23 scales, normalizes them to a common 0–100 range (where 100 = most hydrophobic), and provides tools to:

1. **Define the hydrophobicity space** — visualize how each amino acid's hydrophobicity varies across all scales
2. **Analyze any protein sequence** — plot the per-residue hydrophobicity profile using all 23 scales simultaneously

## Hydrophobicity Distribution

The file `hydrophobicity_distribution.html` is an interactive boxplot showing the normalized hydrophobicity score distribution for each amino acid across all 23 scales, ordered by median hydrophobicity.

## Scripts

### `Scrape_hydrophobicity_data.py`
Scrapes all 23 hydrophobicity scales from ExPASy ProtScale and saves them to `hydrophobicity_scales.csv`.

```bash
python Scrape_hydrophobicity_data.py
```

### `Define_Hydrophobicity_space.py`
Normalizes the scraped scales to a 0–100 range (auto-detecting scale direction) and generates:
- `hydrophobicity_scales_normalized.csv` — normalized scale values
- `hydrophobicity_distribution.html` — interactive boxplot of the hydrophobicity space

```bash
python Define_Hydrophobicity_space.py
```

### `Analyze_Hydrophobicity.py`
Analyzes protein sequences from a FASTA file, plotting per-residue hydrophobicity profiles averaged across all 23 normalized scales. Supports highlighting of codon-level features via an optional highlighted FASTA file with uppercase/lowercase coding.

```bash
python Analyze_Hydrophobicity.py sequences.fasta highlighted_sequences.fasta
```

## Hydrophobicity Scales

The 23 scales included are:

| Scale | Author(s) |
|-------|-----------|
| Hphob.Eisenberg | Eisenberg et al. |
| Hphob.Sweet | Sweet & Eisenberg (OMH) |
| Hphob.Woods | Hopp & Woods |
| Hphob.Manavalan | Manavalan et al. |
| Hphob.Leo | Abraham & Leo |
| Hphob.Black | Black |
| Hphob.Breese | Bull & Breese |
| Hphob.Fauchere | Fauchere et al. |
| Hphob.Guy | Guy |
| Hphob.Janin | Janin |
| Hphob.Miyazawa | Miyazawa et al. |
| Hphob.Argos | Rao & Argos |
| Hphob.Roseman | Roseman |
| Hphob.Tanford | Tanford |
| Hphob.Wolfenden | Wolfenden et al. |
| Hphob.Welling | Welling et al. |
| Hphob.Wilson | Wilson et al. (HPLC) |
| Hphob.Parker | Parker et al. (HPLC) |
| Hphob.pH3.4 | Cowan (HPLC pH 3.4) |
| Hphob.pH7.5 | Cowan (HPLC pH 7.5) |
| Hphob.mobility | Rf mobility |
| Hphob.Chothia | Chothia |
| Hphob.Rose | Rose et al. |

## Installation

```bash
pip install -r requirements.txt
```
