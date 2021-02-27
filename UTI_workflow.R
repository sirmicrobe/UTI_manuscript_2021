#Data analysis for "Performance of conventional urine culture compared to 16S rRNA gene amplicon sequencing in children with suspected urinary tract infection"
#Chris Marshall, PhD; Marcia Kurs-Lasky, MS; Christi L. McElheny, MS; Hui Liu; Nader Shaikh, MD 
#R version 3.6.3
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.10")
library("dada2"); packageVersion("dada2") #1.14.1

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install("phyloseq")
library('phyloseq'); packageVersion('phyloseq') #1.30.0

library(DECIPHER); packageVersion("DECIPHER") # ‘2.8.1’

library("ape"); packageVersion("ape") #5.3 
library("vegan"); packageVersion("vegan") #2.5.6
library("tidyr"); packageVersion("tidyr") #1.1.2


#demultiplexing on beagle with idemp
# -b code = barcode file - Barcode Sampleid; -m n = number of mismatches, default=1; -o output folder
source(idemp -b /home/cwm47/nader_sinusitus/urine_9-2018/dmux_nader_9-2018.txt -I1 /home/cwm47/nader_sinusitus/urine_9-2018/Undetermined_S0_L001_I1_001.fastq -R1 /home/cwm47/nader_sinusitus/urine_9-2018/Undetermined_S0_L001_R1_001.fastq -R2 /home/cwm47/nader_sinusitus/urine_9-2018/Undetermined_S0_L001_R2_001.fastq -m 1 -o /home/cwm47/nader_sinusitus/urine_9-2018/dmux)
#old sinus data
source(idemp -b /home/cwm47/nader_sinusitus/sinus_3-2017/dmux_sinus_3-2017.txt -I1 /home/cwm47/nader_sinusitus/sinus_3-2017/Undetermined_S0_L001_I1_001.fastq.gz -R1 /home/cwm47/nader_sinusitus/sinus_3-2017/Undetermined_S0_L001_R1_001.fastq.gz -R2 /home/cwm47/nader_sinusitus/sinus_3-2017/Undetermined_S0_L001_R2_001.fastq.gz -m 1 -o /home/cwm47/nader_sinusitus/sinus_3-2017/dmux_idemp)


path <- "/home/cwm47/nader_sinusitus/urine_9-2018/dmux"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "*001.fastq_"), `[`, 2)
sample.names[744:851] <- sapply(strsplit(basename(fnFs[744:851]), "*001.fastq.gz_"), `[`, 2)

save.image(file="/home/cwm47/nader_sinusitus/urine_9-2018/dada2/nader_dada2.RData")

#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])


###### Filter and trim ######
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
# filtFs ex. "/home/cwm47/nader_sinusitus/urine_9-2018/dmux/filtered/1_F_filt.fastq.gz" 
# sample.names ex. "1"    "10"   "100"  "1007" "1009" "1010"

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, # maxEE is where the business happens,change strictness here
                     compress=TRUE, multithread=TRUE) #truncLen=c(240,160) use to chop off ends (R1,R2) 
head(out)

#removed samples that had no sequences
filtFs <- filtFs[-c(38,43,57,61,95,196,198,189)]
filtRs <- filtRs[-c(38,43,57,61,95,196,198,189)]
sample.names <- sample.names[-c(38,43,57,61,95,196,198,189)]
which(sample.names=="unsigned") #multiple instances of this sample name so must delete one
sample.names <- sample.names[-822] #

# Learn the error rates - parametric error model
system.time(errF <- learnErrors(filtFs, multithread=TRUE))
system.time(errR <- learnErrors(filtRs, multithread=TRUE))

plotErrors(errF, nominalQ=TRUE) # visualize error rates - black line should follow dots and should decrease with increase quality:  If you are working with a large dataset and the plotted error model does not look like a good fit, you can try increasing the nbases parameter to see if the fit improves.

# Dereplication - combine identical reads into unique seq with abundance value. can be memory intensive. 
derepFs <- derepFastq(filtFs, verbose=TRUE)
system.time(derepRs <- derepFastq(filtRs, verbose=TRUE))

names(derepFs) <- sample.names
names(derepRs) <- sample.names

## Sample Inference - core algorithm of DADA2
dadaFs <- dada(derepFs, err=errF, multithread=TRUE) 
system.time(dadaRs <- dada(derepRs, err=errR, multithread=TRUE)) #this step takes 1492s comp, 377s real time

dadaFs[[1]]
#dada-class: object describing DADA2 denoising results 2 sequence variants were inferred from 23 input unique sequences. Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

#Merge paired end reads - merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region.
system.time(mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)) #194.399s
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
#Error in outs[, 1] : subscript out of bounds -- not sure what this error means: had to upgrade to v1.9.2 here to get around this error

# Construct sequence table of ASVs
system.time(seqtab <- makeSequenceTable(mergers)) #25s
dim(seqtab) #842 samples 11205 ASVs
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab))) # 251nt=58seqs; 252=1442; 253=8481; 254=309
#only use sequences with 240-260 nt
seqtab_orig <- seqtab
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% seq(240,260)]

# Remove chimeras
system.time(seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)) #4s
dim(seqtab.nochim) #842 samples 7930 ASVs; 3275 ASVs removed as chimeric
sum(seqtab.nochim)/sum(seqtab) # frequency of chimeric sequences < 2% (98.3% non chimeric)

# how many reads made it through - tracking reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track,20) #sanity check - make sure you didnt lose a majority of reads at any step

install.packages("seqinr")
library("seqinr"); packageVersion("seqinr") #3.6.1
write.fasta(seqtab.nochim,file.out="seqs_dada2_nochim")

###### Assign taxonomy ####

# naive Bayesian classifier with dada2
#taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v128_train_set.fa.gz", multithread=TRUE) 
# add species level assignments
#taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v128.fa.gz")

# use DECIPHER to classify
library(DECIPHER); packageVersion("DECIPHER") # ‘2.8.1’
#Download the SILVA SSU r132 (modified) file
system(wget http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r132_March2018.RData)

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("/home/cwm47/database/decipher_training_sets/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
system.time(ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)) # use all processors (takes a while)
#use strand="both" if reads not assigned taxonomies (everything is NA NA NA NA)
#user    system   elapsed 
#15902.051   160.796  1993.014 
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
taxa <- taxid
#view taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

save.image(file="/home/cwm47/nader_sinusitus/urine_9-2018/dada2/nader_dada2.RData")

# Evaluate accuracy
unqs.commstk1 <- seqtab.nochim["COMMSTK1",]
unqs.commstk1 <- sort(unqs.commstk1[unqs.commstk1>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.commstk1), "sample sequences present in the Mock community.\n")
#DADA2 inferred 6 sample sequences present in the Mock community.

unqs.pa14 <- seqtab.nochim["Pseudomonoas.aeruginosa.PA14.1.17.17",]
unqs.pa14 <- sort(unqs.pa14[unqs.pa14>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.pa14), "sample sequences present in the Mock community.\n")
#DADA2 inferred 9 sample sequences present in the Mock community.

unqs.s1 <- seqtab.nochim["S1",]
unqs.s1 <- sort(unqs.s1[unqs.s1>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.s1), "sample sequences present in the Mock community.\n")

######################################
#        #### Phyloseq ####          #
######################################

#Import Metadata
samdf_all <- read.csv("sample_names.csv")
rownames(samdf_all) <- samdf_all$sample_names


samdf_urine <- read.csv("../BioSTARRs_the148IDs_sequence_2020.03.11.csv") #updated metadata 3-2020

#create phyloseq object
ps_all <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf_all), 
                   tax_table(taxa))
#otu_table()   OTU Table:         [ 7930 taxa and 842 samples ]

seqspersample <- data.frame(names(sample_sums(ps_all)),sample_sums(ps_all))
write.csv(merge(seqspersample,samdf_all,by="row.names"),"seqs_per_sample.csv")


#subset just the control samples
ps_control <- prune_samples(sample_data(ps_all)$project == "control", ps_all) # Remove mock sample
ps_control #otu_table()   OTU Table:         [ 7930 taxa and 14 samples ]
sum(taxa_sums(ps_control) == 0)
ps_control = prune_taxa(taxa_sums(ps_control) > 0, ps_control) #[ 188 taxa and 14 samples ]

#pseudomonas only
control_pa14 <- prune_samples(sample_data(ps_control)$sample_names == "Pseudomonoas.aeruginosa.PA14.1.17.17", ps_control)
control_pa14 = prune_taxa(taxa_sums(control_pa14) > 0, control_pa14) #[ 9 taxa and 1 samples ]

ps_urine <- prune_samples(sample_data(ps_all)$project == "urine", ps_all) # Remove mock sample
ps_urine #otu_table()   OTU Table:         [ 7930 taxa and 134 samples ]
ps_urine_orig <- ps_urine
#prune samples with less than 1000 reads
ps_urine = prune_samples(sample_sums(ps_urine)>=1000, ps_urine)
sum(taxa_sums(ps_urine) == 0)
ps_urine = prune_taxa(taxa_sums(ps_urine) > 0, ps_urine) #[ 864 taxa and 122 samples ]

#create phylogenetic tree
random_tree_urine = rtree(ntaxa(ps_urine), rooted=TRUE, tip.label=taxa_names(ps_urine))
plot(random_tree_urine)
#merge into one phyloseq object
samdf_urine <- merge(data.frame(sample_data(ps_urine)), samdf_urine,by="sample_names")
row.names(samdf_urine) <- as.character(samdf_urine$sample_names)
sam_urine <- sample_data(samdf_urine)

sample_data(ps_urine) <- sam_urine
sample_data(ps_urine)$sample_names <- as.character(sample_data(ps_urine)$sample_names)
sample_data(ps_urine)$growth_paper <- as.factor(as.character(toupper(sample_data(ps_urine)$growth_paper)))
#drop samples 6454, 9514, 9522, 9523
ps_urine <- subset_samples(ps_urine, !(sample_names %in% c("6454", "9514", "9522","9523"))) #118 samples

#output tax table
getwd()
TT_urine <- as(tax_table(ps_urine),"matrix")
head(TT_urine)
write.csv(as.data.frame(TT_urine),file="urine_tax_table.csv",quote=F)
#output otu_table
write.csv(otu_table(ps_urine),"urine_otu_table.csv",quote=F)


#### barplots ####
ps_urine.prop <- transform_sample_counts(ps_urine, function(OTU) OTU/sum(OTU))

plot_bar(ps_urine.prop , x="sample_names", fill="family")+ scale_fill_manual(values = bar_colors) +coord_flip()+facet_wrap(~uti_dx,scales="free_y")#+theme(legend.position = "none")

#top30
top30 <- names(sort(taxa_sums(ps_urine), decreasing=TRUE))[1:30]
ps_urine.top30 <- prune_taxa(top30, ps_urine.prop)

# Figure 1 barplots
plot_bar(ps_urine.top30 , x="sample_names", fill="family")+ scale_fill_manual(values = c("grey26","red","blue4","royalblue","mediumpurple1","darkorange","cyan2", "darkgreen","chartreuse3", "mediumorchid3","orangered3","salmon1","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "yellow1", "#599861","brown2","darkcyan","orchid1","orange1","skyblue4","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4")) +facet_wrap(~uti_dx,scales="free_y")+coord_flip()

#### plot alpha diversity ####

# Figure 2A = Shannon diversity based on clinical score
p_rich_urine3 <- plot_richness(ps_urine,x="uti_dx", measures=c("Shannon"))
p_rich_urine3$layers <- p_rich_urine3$layers[-1]
p_rich_urine3 + geom_boxplot() + geom_point()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
# Figure 2B = Shannon diversity based on urine collection method
p_rich_urine4 <- plot_richness(ps_urine,x="Method.of.urine.collection", measures=c("Shannon"))
p_rich_urine4$layers <- p_rich_urine4$layers[-1]
p_rich_urine4 + geom_boxplot() + geom_point()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#stats
#generate richness data frame
d_rich <- estimate_richness(ps_urine, measures = "Shannon")
d_rich$clean_catch <- sample_data(ps_urine)$Method.of.urine.collection
d_rich$uti_dx <- sample_data(ps_urine)$uti_dx

#Kruskal wallis test
kruskal.test(Shannon ~ uti_dx, data = d_rich) #Kruskal-Wallis chi-squared = 55.757, df = 5, p-value = 9.12e-11
library("FSA") #v0.8.30 Ogle, D.H., P. Wheeler, and A. Dinno. 2020. FSA: Fisheries Stock Analysis. R package version 0.8.30,https://github.com/droglenc/FSA.
dunnTest(Shannon ~ as.factor(uti_dx),
         data=d_rich,
         method="bh") #multiple comparison test 

#independent 2-group Mann-Whitney U Test
wilcox.test(Shannon ~ clean_catch, data = d_rich) #W = 1005.5, p-value = 8.347e-05
t.test(Shannon ~ clean_catch, data = d_rich) #t = -3.9958, df = 96.141, p-value = 0.0001263



#### beta diversity NMDS ordination ####

#Figure S2
#Transform data to proportions as appropriate for Bray-Curtis distances
ps_urine.prop <- transform_sample_counts(ps_urine, function(otu) otu/sum(otu))
#Bray-Curtis dissimilarity for NMDS ordination
ps_urine.ord.nmds.bray <- ordinate(ps_urine.prop, method="NMDS", distance="bray")
#plot ordination
plot_ordination(ps_urine.prop, ps_urine.ord.nmds.bray, color="uti_dx",title="Bray NMDS") +geom_point(size=4) +scale_color_manual(values=c("#9750a1","#6fac5d","#bc7d39")) 

#testing for significance:
d = distance(ps_urine.prop, "bray")
myDataadonis = adonis2(d ~ uti_dx, as(sample_data(ps_urine.prop), "data.frame"))
myDataadonis 
#looking for post hoc from adonis
library("RVAideMemoire") #v 0.9-78
pairwise.perm.manova(d,as(sample_data(ps_urine.prop), "data.frame")$uti_dx,nperm=999)
