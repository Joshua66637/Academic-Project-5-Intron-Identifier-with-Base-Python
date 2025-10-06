# Split Read Junction Extractor

## Overview
This project identifies **splice junctions** from RNA-seq alignments in **SAM** format and associates them with genes using a provided annotation table.  
Split reads (those with `N` in their CIGAR strings) indicate intron locations. The script counts how many reads support each junction and outputs all unique junctions per gene.

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

## How the Script Works

1. Skips SAM header lines (`@`).  
2. Keeps only uniquely aligned reads (`NH:i:1`).  
3. Detects split reads using `N` in the CIGAR string.  
4. Calculates each intron’s start and end positions.  
5. Groups identical junctions and counts supporting reads.  
6. Matches junctions to genes and writes results to `XXXXXXXX.txt` (replace with your student number).  
   Each gene’s junctions are followed by a blank line.

---

## Usage

```
python3 myScript.py mySamFile.sam myInputTable.txt
```
Example:
```
python3 myScript.py alignments.sam genes.txt
```
