# Split Read Junction Extractor

## Overview
This project identifies **splice junctions** from RNA-seq alignments in **SAM** format and associates them with genes using a provided annotation table.  
Split reads are caused due to introns in the Genomic Reference Sequence which are absent in the mRNA. These read span junctions and have `N`s in their CIGAR strings. These indicate intron locations. The script counts how many reads support each junction and outputs all unique junctions per gene.

---

## Input Files

1. **SAM file** – RNA-seq alignments.  
   Only reads with `NH:i:1` (unique alignments) are processed.  
   Required fields: chromosome, start position, CIGAR string, and the final `NH:i:x` tag.

2. **Gene table file** – Tab-separated text file with:  
   - Gene ID  
   - Transcript ID  
   - Genomic location (`chromosome:start..end(strand)`)

---

## Usage

```
python3 IntronFinder.py Samfile.sam GeneLocationTable.txt
```
