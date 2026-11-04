# 🧬 Introduction to Genomic Analysis with Bioconductor (R Workshop)

## 📌 About This Workshop

This workshop was organized by **Nikolaos Katsikavalis** and the **IEEE Student Branch Lamia** at the **Department of Computer Science and Biomedical Informatics, University of Thessaly**.

This repository contains the complete hands-on code, documentation, and practical workflow from the event. The primary objective was to introduce students and researchers to genomic data analysis in R using the **Bioconductor** ecosystem—covering everything from core `S4` object architecture and sequence manipulation to genomic ranges (`GRanges`), transcriptomic annotations (`TxDb`), and automated Next-Generation Sequencing (NGS) quality control (`Rqc`).

* **Speaker:** Nikolaos Katsikavalis
* **Organized by:** Nikolaos Katsikavalis & IEEE Student Branch Lamia
* **Affiliation:** Department of Computer Science and Biomedical Informatics, University of Thessaly
* **Contact:** [nkatsikavalis@uth.gr](mailto:nkatsikavalis@uth.gr) 
---

## 📂 Repository Contents

* `Bioconductor_Workshop_Script.R`: The complete, self-contained R script with all workshop modules (0 to 9).
* `Οδηγός_Workshop_Bioconductor.pdf`: Comprehensive step-by-step workshop guide and documentation (in Greek).
* `image_f1e526.png`: Per-sample FASTQ mean quality Boxplot generated during the QC module.

---

## 🚀 Quick Start & Automated Setup (Module 0)

To ensure 100% reproducibility and fast execution during the live session without downloading massive multi-gigabyte NGS files, the script automatically generates a local `data/` directory containing lightweight mock **FASTA** (`mock_virus.fa`) and **FASTQ** (`sample_1.fq`, `sample_2.fq`) files, and installs all required Bioconductor packages.

```R
if (!dir.exists("data")) dir.create("data")

# Generate mock FASTA & FASTQ files
fasta_content <- c(">seq1", "ATGCGTACGTAGCTAG", ">seq2", "CGTACGTAGCTAGCTA")
writeLines(fasta_content, "data/mock_virus.fa")

fastq1_content <- c("@SEQ1", "GATTTGGGGTTCAAAG", "+", "!''*((((***+))%%", 
                    "@SEQ2", "ATCGTACGATCGATCG", "+", "IIIIIIIIIIIIIIII")
writeLines(fastq1_content, "data/sample_1.fq")

fastq2_content <- c("@SEQ1_S2", "GATTTGGGGTTCAAAG", "+", "!''*((((***+))%%", 
                    "@SEQ2_S2", "ATCGTACGATCGATCG", "+", "IIIIIIIIIIII####")
writeLines(fastq2_content, "data/sample_2.fq")

# Install required Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "Biostrings", 
  "GenomicRanges", 
  "Rqc", 
  "BSgenome.Scerevisiae.UCSC.sacCer3",
  "TxDb.Hsapiens.UCSC.hg38.knownGene"
), update = FALSE, ask = FALSE)
```

---

## 🔬 Workshop Modules & Expected Outputs

### 1. Environment & Validity Check (`BiocManager`)
Verifies that the R environment and installed Bioconductor packages are properly synchronized.
```R
library(BiocManager)
BiocManager::version()
BiocManager::valid()
sessionInfo()
```

### 2. S4 Architecture & `S4Vectors` (`DataFrame`)
Explores Bioconductor's strict S4 class system using `DataFrame` (capital D), which safely stores complex genomic structures in predefined slots.
```R
library(S4Vectors)
df <- DataFrame(chromosome = c("chr1", "chr2"), genes = c(10, 15))
isS4(df)
showMethods(class = "DataFrame")
```
* **Expected Output:** `isS4(df)` returns `[1] TRUE`.

### 3. Biological Sequence Manipulation (`Biostrings`)
Handles DNA/RNA sequences, translates codons (triplets) into amino acids (`*` indicates a Stop codon), computes the 5' -> 3' reverse complement, and parses FASTA files.
```R
library(Biostrings)
dna_seq <- DNAString("ATGCGTACGTAGCTAGCTAGT")
rna_seq <- RNAString(dna_seq)
protein_seq <- translate(dna_seq)
reverseComplement(dna_seq)

fasta_reads <- readDNAStringSet("data/mock_virus.fa")
width(fasta_reads)
```
* **Expected Output:**
  ```text
  > rna_seq
  21-letter RNAString object
  seq: AUGCGUACGUAGCUAGCUAGU
  > protein_seq
  7-letter AAString object
  seq: MRT*LAS
  > reverseComplement(dna_seq)
  21-letter DNAString object
  seq: ACTAGCTAGCTACGTACGCAT
  > width(fasta_reads)
  [1] 16 16
  ```

### 4. Reference Genomes (`BSgenome`)
Loads the full reference genome of *Saccharomyces cerevisiae* (yeast), inspects chromosome lengths, and extracts specific genomic coordinates from the mitochondrial chromosome (`chrM`).
```R
library(BSgenome.Scerevisiae.UCSC.sacCer3)
yeast <- BSgenome.Scerevisiae.UCSC.sacCer3
seqlengths(yeast)
chrM_sub <- getSeq(yeast, names = "chrM", start = 1, end = 50)
chrM_sub
```
* **Expected Output:**
  ```text
  50-letter DNAString object
  seq: TTCATAATTAATTTTTTATATATATATTATATTATAATATTAATTTATAT
  ```

### 5. Genomic Intervals (`IRanges` & `GRanges`)
Demonstrates the difference between simple 1-based closed integer intervals (`IRanges`) and biological genomic coordinates (`GRanges`) that include chromosome (`seqnames`) and strand direction (`+` forward 5' -> 3', `-` reverse 3' -> 5').
```R
library(IRanges)
library(GenomicRanges)

ir <- IRanges(start = c(1, 100), end = c(50, 250))
gr <- GRanges(
  seqnames = c("chr1", "chr2"),
  ranges = IRanges(start = c(100, 500), end = c(300, 750)),
  strand = c("+", "-")
)
gr
```
* **Expected Output:**
  ```text
  GRanges object with 2 ranges and 0 metadata columns:
        seqnames    ranges strand
           <Rle> <IRanges>  <Rle>
    [1]     chr1   100-300      +
    [2]     chr2   500-750      -
  ```

### 6. Spatial Interval Overlaps (`findOverlaps`)
Identifies spatial intersections between experimental regions (`query`) and reference regions (`subject`) using the 1-based closed interval formula (End = Start + Width - 1).
```R
query <- GRanges("chr1", IRanges(start = c(50, 150), width = 60))     # 50-109, 150-209
subject <- GRanges("chr1", IRanges(start = c(100, 200), width = 80))  # 100-179, 200-279

hits <- findOverlaps(query, subject)
hits
subsetByOverlaps(query, subject)
```
* **Expected Output:** Returns `3 hits` (`Query 1 -> Subject 1`, `Query 2 -> Subject 1`, and `Query 2 -> Subject 2`).

### 7. Genome Annotations & Transcriptomes (`TxDb`)
Queries the Human Genome annotation database (`UCSC hg38`), restricts the scope to chromosome X (`chrX`), and groups exons by transcript (`tx`).
```R
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(txdb) <- "chrX"
exons_per_tx <- exonsBy(txdb, by = "tx")
exons_per_tx[1:2]
```
* **Expected Output:** Returns a `GRangesList` containing the exon coordinates for transcripts `372516` (2 exons) and `372517` (8 exons) on `chrX`.

### 8. Reading & Filtering NGS Reads (`ShortRead`)
Parses raw FASTQ sequencing files, inspects DNA sequences (`sread`) and ASCII-encoded Phred quality scores (`quality`), and applies a custom minimum length filter (>= 16 sequencing cycles/bases).
```R
library(ShortRead)
fq <- readFastq("data/sample_1.fq")
sread(fq)
quality(fq)

len_filter <- srFilter(function(x) width(x) >= 16, name = "MinLength16")
fq_filtered <- fq[len_filter(fq)]
length(fq_filtered)
```
* **Expected Output:**
  ```text
  class: ShortReadQ
  length: 2 reads; width: 16 cycles
  [1] 2
  ```

### 9. Automated FASTQ Quality Control (`Rqc`)
Performs batch quality control across all `.fq` files in the `data/` directory and visualizes the per-sample mean quality distribution using a Boxplot.
```R
library(Rqc)
qc_results <- rqc(path = "data", pattern = ".fq$", sample = TRUE)
rqcReadQualityBoxPlot(qc_results)
```

#### 📊 Quality Control Visualization (Boxplot)
![Rqc Read Quality Boxplot](Rplot.png)

* **`sample_1.fq`:** Exhibits higher overall read quality (median around **Q24**, reaching up to **Q40** for `I` bases), with wider variance due to lower quality scores in the initial bases of the first read (`!''*`).
* **`sample_2.fq`:** Shows a lower quality distribution (median around **Q19**), directly reflecting the sharp quality drop at the 3' tail of the second read (`####`).
