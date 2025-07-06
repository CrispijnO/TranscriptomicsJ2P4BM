#verander de huidige werkmap naar de locatie met de data erin
setwd("/data/locatie")
#installatie benodigde packages
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
BiocManager::install("scales")
#inladen van deze packages
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
library(scales)

#indexeren van het menselijk genoom, waar basename de bestand naam bepaald, reference de referentie is waar alle genen tot behoren, memory bepaalt hoeveel RAM het proces mag gebruiken. indexsplit splitst de referentie en index om de functie te versnellen
buildindex(
  basename = "ref_human",
  reference = "Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
  memory = 16000,
  indexSplit = TRUE
)
#alignen van de sample bestanden op de gebouwde index, index refereert de index die wordt gebruikt, readfile specificeert welk bestand wordt gebruikt om te alignen en output file is de bestandnaam van het aligned sample
alignCtrl1 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785819_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785819_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785819_control.BAM")
alignCtrl2 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785820_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785820_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785820_control.BAM")
alignCtrl3 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785828_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785828_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785828_control.BAM")
alignCtrl4 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785831_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785831_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785831_control.BAM")
alignCase1 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785979_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785979_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785979_case.BAM")
alignCase2 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785980_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785980_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785980_case.BAM")
alignCase3 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785986_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785986_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785986_case.BAM")
alignCase4 = align(index = "ref_human", readfile1 = "Data_RA_raw/SRR4785988_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785988_2_subset40k.fastq", output_file = "Data_RA_raw/aligned/SRR4785988_case.BAM")

#zet de bestand namen de aligned samples bestand in een variabele
allSamples = c("Data_RA_raw/aligned/SRR4785819_control.BAM", "Data_RA_raw/aligned/SRR4785820_control.BAM", "Data_RA_raw/aligned/SRR4785828_control.BAM", "Data_RA_raw/aligned/SRR4785831_control.BAM", "Data_RA_raw/aligned/SRR4785979_case.BAM", "Data_RA_raw/aligned/SRR4785980_case.BAM", "Data_RA_raw/aligned/SRR4785986_case.BAM", "Data_RA_raw/aligned/SRR4785988_case.BAM")

#maakt een count matrix van de sample bestanden, files specifeert welke bestanden gebruikt worden om te tellen, annot.ext is het annotatie bestand bestaande uit gensequenties en gen naam, dat gebruikt wordt om te tellen hoe vaak een gen voorkomt.isGTFAnnotationFile specifeert het annotatie bestand format. isPairedEnd TRUE is nodig als het te tellen bestand dubbelstrengs is. GTF.attrType specifeert de gen naam/id dat als naam gebruikt moet worden als het gen geteld is. useMetaFeatures noemt de rijnamen van de getelde gen naar het gen id
count_matrix = featureCounts(
  files = allSamples,
  annot.ext = "Homo_sapiens.GRCh38.114.gtf.gz",
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

#leest het count_matrix.txt bestand dat nodig was en zet dit in een variabele
counts = read.table("count_matrix.txt")
#verandert de naam van de kolommen naar het sample id
colnames(counts) = c("SRR4785819_control", "SRR4785820_control", "SRR4785828_control", "SRR4785831_control", "SRR4785979_case", "SRR4785980_case", "SRR4785986_case", "SRR4785988_case")
# maakt een data frame waarin de status van het elk monster wordt toegewezen zodat dit gebruikt kan worden in de DESeqDataSetFromMatrix functie
treatment = c("control", "control", "control", "control", "case", "case", "case", "case")
treatment_table = data.frame(treatment)
rownames(treatment_table) = c("SRR4785819_control", "SRR4785820_control", "SRR4785828_control", "SRR4785831_control", "SRR4785979_case", "SRR4785980_case", "SRR4785986_case", "SRR4785988_case")

#Zet de count_matrix tabel om in een DESeq data set met juiste sample namen en de RA of controle status
dds = DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = treatment_table,
  design = ~ treatment
)
#voert de differentiële genexpressie analyse uit
dds = DESeq(dds)
#zet de resultaten van de analyse in een tabel en in een variabele waar de case case wordt vergeleken tegen de controle 
ddsResults = results(dds, contrast = c("treatment", "case", "control"))
#haalt lege resultaten uit de tabel
ddsResults = na.omit(ddsResults)
#zet de differentiële genexpressie analyse resultaten tabel in een apart bestand, waar eerst de tabel wordt gespecificeerd. file bepaald de bestandsnaam van het nieuwe bestand. row.names specificeert dat de tabel rij namen heeft, hetzelfde voor col names. sep specificeert de seperator tussen elke gegevens in de tabel
write.table(ddsResults, file = "Differentiële genexpressie analyse resulaten.csv", row.names = TRUE, col.names = TRUE, sep= ";")
#maakt de volcano plot uit ddsresults. lab labeld de namen van de data punten, x specificeert de naam van de x-as, hetzelfde voor y. xlim bepaalt het de uiterste waarde van de x-as
EnhancedVolcano(
  ddsResults,
  lab = rownames(ddsResults),
  x = "log2FoldChange",
  y = "padj",
  xlim = c(-13,13)
)
#filtert de differentiële genexpressie analyse op een p van <0.05 en een fold change van >1 of <-1
ddsResultsF = ddsResults[ddsResults$padj < 0.05 & (ddsResults$log2FoldChange < -1 | ddsResults$log2FoldChange > 1),]
#zet de namen van de gefilterde genen in een variabele
genesToTest = rownames(ddsResultsF)
#voert de GO enrichment analyse uit met de genesToTest variabele. OrgDb bepaalt de annotatie die wordt gebruikt. Hier is het de annotatie van de mens. keyType specificeert de format van de gennamen, ont bepaalt het ontologie niveau
GOResBP = enrichGO(genesToTest, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
#zet de GO enrichment resultaten in een apart bestand, waar eerst de tabel wordt gespecificeerd. file bepaald de bestandsnaam van het nieuwe bestand. row.names specificeert dat de tabel rij namen heeft, hetzelfde voor col names. sep specificeert de seperator tussen elke gegevens in de tabel
write.table(GOResBP@result, file = "GO Enrichment resultaten.csv", row.names = TRUE, col.names = TRUE, sep= ";")
#sorteert de data data op hoog naar laag
GOResBPplot = GOResBP@result[order(1:12),]
#maakt een grafiek van de GO Enrichment analyse. aes(wat er op x-as staat, wat er op de y-as staat) + geom_col specificeert wat voor grafiek het wordt(fil = kleur van de balken) + labs(x= naam x-as, y= naam y-as, title= naam van de grafiek) + theme_minimal veranderd de stijl van de grafiek + guides(fill=none) is nodig voor het veranderen van de letter grootte + scale_y_discrete veranderd de lettergrootte van de y-as(labels = label_wrap is geeft het maximaal aantal karakters voor de labels van de y-as, functie van scales)
ggplot(GOResBPplot, aes(-log10(p.adjust), reorder(Description, -log10(p.adjust)))) + geom_col(fill = "cornflowerblue") + labs(x = "-log10 P Waarde", y = "GO-Term", title = "Meeste verrijkte GO-termen met de hooste significantie") + theme_minimal() + guides(fill = "none") + scale_y_discrete(labels = label_wrap(40))
#verwijderd onnodige resultaten voor de KEGG pathway analyse
ddsResultsF[1] = NULL
ddsResultsF[2:5] = NULL

#voert de KEGG pathway analyse uit. gene.data geeft de data voor analyse, species specificeert om welk organisme het gaat. pathway.id specificeert het te gebruiken pathway met id van https://www.kegg.jp/ . gene.idtype specificeert dat de gennamen in de data in het SYMBOL format staat, limit geeft het minimale tot maximale kleur verandering die de functie aan een gen kan geven
pathview = pathview(
  gene.data = ddsResultsF,
  species = "hsa",
  pathway.id = "hsa05323",
  gene.idtype = "SYMBOL",
  limit = list(gene = 5)
)