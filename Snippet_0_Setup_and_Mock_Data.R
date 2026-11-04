# ==============================================================================
# Bioconductor Workshop: Complete Executable R Script
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Προετοιμασία Περιβάλλοντος, Mock Δεδομένα & Εγκατάσταση Πακέτων
# ------------------------------------------------------------------------------
# Ορισμός του φακέλου εργασίας (αφαιρέστε το σχόλιο αν χρειάζεται)
# setwd("~/Desktop/Bioconductor_Workshop")

if (!dir.exists("data")) dir.create("data")

message("Δημιουργία δεδομένων...")
fasta_content <- c(">seq1", "ATGCGTACGTAGCTAG", ">seq2", "CGTACGTAGCTAGCTA")
writeLines(fasta_content, "data/mock_virus.fa")

fastq1_content <- c("@SEQ1", "GATTTGGGGTTCAAAG", "+", "!''*((((***+))%%", "@SEQ2", "ATCGTACGATCGATCG", "+", "IIIIIIIIIIIIIIII")
writeLines(fastq1_content, "data/sample_1.fq")

fastq2_content <- c("@SEQ1_S2", "GATTTGGGGTTCAAAG", "+", "!''*((((***+))%%", "@SEQ2_S2", "ATCGTACGATCGATCG", "+", "IIIIIIIIIIII####")
writeLines(fastq2_content, "data/sample_2.fq")

message("Εγκατάσταση πακέτων...")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "Biostrings", 
  "GenomicRanges", 
  "Rqc", 
  "BSgenome.Scerevisiae.UCSC.sacCer3",
  "TxDb.Hsapiens.UCSC.hg38.knownGene"
), update = FALSE, ask = FALSE)
message("Όλα έτοιμα!")

# ------------------------------------------------------------------------------
# 1. Έλεγχος Έκδοσης Bioconductor & Περιβάλλοντος
# ------------------------------------------------------------------------------
library(BiocManager)
BiocManager::version()
BiocManager::valid()
sessionInfo()

# ------------------------------------------------------------------------------
# 2. Αρχιτεκτονική S4 & S4Vectors (DataFrame)
# ------------------------------------------------------------------------------
library(S4Vectors)

# Δημιουργία ενός DataFrame (αντικείμενο S4) με δεδομένα
df <- DataFrame(chromosome = c("chr1", "chr2"), genes = c(10, 15))
df

# Έλεγχος αν το αντικείμενο ανήκει στην αρχιτεκτονική S4
isS4(df)

# Εμφάνιση των διαθέσιμων μεθόδων για την κλάση DataFrame
showMethods(class = "DataFrame")

# ------------------------------------------------------------------------------
# 3. Χειρισμός Αλληλουχιών με το Biostrings
# ------------------------------------------------------------------------------
library(Biostrings)

# Ορισμός αλληλουχίας DNA
dna_seq <- DNAString("ATGCGTACGTAGCTAGCTAGT")

# Μεταγραφή σε RNA και μετάφραση σε πρωτεΐνη
rna_seq <- RNAString(dna_seq)
rna_seq

protein_seq <- translate(dna_seq)
protein_seq

# Συμπληρωματικότητα και αντίστροφη φορά (5' -> 3')
reverseComplement(dna_seq)

# Ανάγνωση αρχείου FASTA από τον φάκελο data
fasta_reads <- readDNAStringSet("data/mock_virus.fa")
fasta_reads

width(fasta_reads)

# ------------------------------------------------------------------------------
# 4. Πλήρη Γονιδιώματα Αναφοράς με το BSgenome
# ------------------------------------------------------------------------------
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Ανάθεση του γονιδιώματος σε μεταβλητή
yeast <- BSgenome.Scerevisiae.UCSC.sacCer3

# Επισκόπηση χρωμοσωμάτων και μεγεθών
seqlengths(yeast)

# Εξαγωγή μιτοχονδριακής αλληλουχίας (chrM, βάσεις 1 έως 50)
chrM_sub <- getSeq(yeast, names = "chrM", start = 1, end = 50)
chrM_sub

# ------------------------------------------------------------------------------
# 5. Διαστήματα & Γονιδιωματικές Συντεταγμένες (IRanges & GRanges)
# ------------------------------------------------------------------------------
library(IRanges)
library(GenomicRanges)

# Βασικά διαστήματα μήκους (IRanges)
ir <- IRanges(start = c(1, 100), end = c(50, 250))
ir

# Γονιδιωματικά διαστήματα με χρωμόσωμα και αλυσίδα (GRanges)
gr <- GRanges(
  seqnames = c("chr1", "chr2"),
  ranges = IRanges(start = c(100, 500), end = c(300, 750)),
  strand = c("+", "-")
)
gr

# ------------------------------------------------------------------------------
# 6. Επικαλύψεις Περιοχών (findOverlaps)
# ------------------------------------------------------------------------------
# Ορισμός περιοχών
query <- GRanges("chr1", IRanges(start = c(50, 150), width = 60))
subject <- GRanges("chr1", IRanges(start = c(100, 200), width = 80))

# Εντοπισμός επικαλύψεων
hits <- findOverlaps(query, subject)
hits

# Φιλτράρισμα περιοχών με επικάλυψη
subsetByOverlaps(query, subject)

# ------------------------------------------------------------------------------
# 7. Γονιδιωματικός Σχολιασμός & Μεταγραφώματα (TxDb)
# ------------------------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Περιορισμός της ανάλυσης στο χρωμόσωμα X
seqlevels(txdb) <- "chrX"

# Εξαγωγή εξονίων ανά μετάγραφημα
exons_per_tx <- exonsBy(txdb, by = "tx")
exons_per_tx[1:2]

# ------------------------------------------------------------------------------
# 8. Ανάγνωση & Φιλτράρισμα NGS Reads (ShortRead)
# ------------------------------------------------------------------------------
library(ShortRead)

# Ανάγνωση αρχείου FASTQ
fq <- readFastq("data/sample_1.fq")

# Εξαγωγή αλληλουχιών DNA και σκορ ποιότητας Phred
sread(fq)
quality(fq)

# Δημιουργία φίλτρου ελάχιστου μήκους (16 κύκλοι / βάσεις)
len_filter <- srFilter(function(x) width(x) >= 16, name = "MinLength16")

# Εφαρμογή φίλτρου και έλεγχος επιβίωσης reads
fq_filtered <- fq[len_filter(fq)]
fq_filtered
length(fq_filtered)

# ------------------------------------------------------------------------------
# 9. Ποιοτικός Έλεγχος Αρχείων FASTQ (Rqc)
# ------------------------------------------------------------------------------
library(Rqc)

# Αυτόματη σάρωση του φακέλου data για αρχεία .fq
qc_results <- rqc(path = "data", pattern = ".fq$", sample = TRUE)

# Δημιουργία γραφήματος boxplot συνολικής ποιότητας ανά δείγμα
rqcReadQualityBoxPlot(qc_results)