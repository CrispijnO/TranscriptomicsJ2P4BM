setwd("/Users/crisp/Downloads/School/Transcriptomics/")
install.packages('BiocManager')
BiocManager::install("Rsubread")
BiocManager::install("DESeq2")
BiocManager::install("KEGGREST")
BiocManager::install("EnhancedVolcano")
BiocManager::install("pathview")
BiocManager::install("AnnotationDbi")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ggplot2")
BiocManager::install("dplyr")
BiocManager::install("scales")
library(BiocManager)
library(Rsubread)
library(DESeq2)
library(KEGGREST)
library(EnhancedVolcano)
library(pathview)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(scales)

buildindex(
  basename = "ref_human",
  reference = "Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
  memory = 16000,
  indexSplit = TRUE
)

alignCtrl1 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785819_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785819_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785819_control.BAM")
alignCtrl2 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785820_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785820_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785820_control.BAM")
alignCtrl3 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785828_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785828_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785828_control.BAM")
alignCtrl4 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785831_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785831_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785831_control.BAM")
alignCase1 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785979_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785979_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785979_case.BAM")
alignCase2 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785980_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785980_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785980_case.BAM")
alignCase3 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785986_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785986_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785986_case.BAM")
alignCase4 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785988_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785988_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785988_case.BAM")

allSamples = c("Data_RA_raw/aligned/SRR4785819_control.BAM", "Data_RA_raw/aligned/SRR4785820_control.BAM", "Data_RA_raw/aligned/SRR4785828_control.BAM", "Data_RA_raw/aligned/SRR4785831_control.BAM", "Data_RA_raw/aligned/SRR4785979_case.BAM", "Data_RA_raw/aligned/SRR4785980_case.BAM", "Data_RA_raw/aligned/SRR4785986_case.BAM", "Data_RA_raw/aligned/SRR4785988_case.BAM")

count_matrix = featureCounts(
  files = allSamples,
  annot.ext = "Homo_sapiens.GRCh38.114.gtf.gz",
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

counts = read.table("count_matrix.txt")
colnames(counts) = c("SRR4785819_control", "SRR4785820_control", "SRR4785828_control", "SRR4785831_control", "SRR4785979_case", "SRR4785980_case", "SRR4785986_case", "SRR4785988_case")

treatment = c("control", "control", "control", "control", "case", "case", "case", "case")
treatment_table = data.frame(treatment)
rownames(treatment_table) = c("SRR4785819_control", "SRR4785820_control", "SRR4785828_control", "SRR4785831_control", "SRR4785979_case", "SRR4785980_case", "SRR4785986_case", "SRR4785988_case")

dds = DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = treatment_table,
  design = ~ treatment
)
dds = DESeq(dds)
ddsResults = results(dds, contrast = c("treatment", "case", "control"))
ddsResults = na.omit(ddsResults)
write.table(ddsResults, file = "Differentiële genexpressie analyse resulaten.csv", row.names = TRUE, col.names = TRUE, sep= ";")

EnhancedVolcano(
  ddsResults,
  lab = rownames(ddsResults),
  x = "log2FoldChange",
  y = "padj",
  xlim = c(-13,13)
)

ddsResultsF = ddsResults[ddsResults$padj < 0.05 & (ddsResults$log2FoldChange < -1 | ddsResults$log2FoldChange > 1),]
genesToTest = rownames(ddsResultsF)
GOResBP = enrichGO(genesToTest, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
write.table(GOResBP@result, file = "GO Enrichment resultaten.csv", row.names = TRUE, col.names = TRUE, sep= ";")
GOResBPplot = GOResBP@result[order(1:12),]

ggplot(GOResBPplot, aes(-log10(p.adjust), reorder(Description, -log10(p.adjust)))) + geom_col(fill = "cornflowerblue") + labs(x = "-log10 P Waarde", y = "GO-Term", title = "Meeste verrijkte GO-termen met de hooste significantie") + theme_minimal() + guides(fill = "none") + scale_y_discrete(labels = label_wrap(40))

ddsResultsF[1] = NULL
ddsResultsF[2:5] = NULL

pathview = pathview(
  gene.data = ddsResultsF,
  species = "hsa",
  pathway.id = "hsa05323",
  gene.idtype = "SYMBOL",
  limit = list(gene = 5)
)