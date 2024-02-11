##Sarracenia purpurea Phytotelmata Metagenomics Project
##Paul Melchior, Depts of Biology, NHCC/Bemidji State University; Dept of Vet. Med., Univ. of Minnesota 2024

#The following R code was used to analyze 16S Illumina 2 x 150 amplicons from 94 samples collected from multiple plants in 
#three Minnesota bogs during one growing season.

##############################################################################################################################3

#PART I:  Filter, trim, and infer DNA sequences from 16S Illumina amplicons using Cutadapt and dada2 packages in R.

#I.1  Removal of primers with Cutadapt (Python)

#Preparation of Files:
  
ls *_R1_001.fastq.gz

#Generate a (file) count of all *.gz files
ls *_R1_001.fastq.gz | wc -l

#Create a file with the prefix names of every fastq sample with the R1/R2 component removed.  Use the python sed #command and 's' operator to shorten files names (remove _R1_001.fastq.gz, replace with nothin')
#Then use a pipe | to output these file names to a new file called 'Sample_names' 

ls *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' > Sample_names
cat Sample_names

#Running Cutadapt in Conda:

#Activate conda environment
conda activate

#Show conda environments available
conda info --envs

#Activate cutadapt environment
conda activate cutadaptenv

#Run Cutadapt without min and max length info.  See Astrobiomike page for that
#Code includes the forward and reverse primer sequences AND their reverse compliments
#Note:  GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC are the R1 forward and rev comp of reverse 16s primer sequences. #The other two are the R2 (thus 'reverse' and rev comp of forward 16s primers
for sample in $(cat Samples_names); do echo "On sample: $sample"; cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC -A ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC -o ${sample}_R1_001_trimmed.fastq.gz -p ${sample}_R2_001_trimmed.fastq.gz ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz -m 1 >> cutadapt_primer_trimming_stats.txt 2>&1; done

#Output files will have sample name, run, etc, and *.trimmed.fastqz

##################################################################

#I.2  DADA2:  Filter, trim, and infer DNA sequences with dada2

library(dada2)

#a. Set path to directory with cutadapt-trimmed fastq files

path <-"C:/Users/sj4879dj/OneDrive - MNSCU/MelchiorFiles/Research/Research - Sarrecenia Project/Sarr_2023/SarrMN_Final/CAtrimmed" 

#b. Create matched list of file names for manipulation:
  # Forward and reverse fastq filenames have format: SAMPLENAME_R1_001_trimmed.fastq 
  # and SAMPLENAME_R2_001_trimmed.fastq

fnFs <- sort(list.files(path, pattern="_R1_001_trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001_trimmed.fastq", full.names = TRUE))

#c. Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

#d. Plot and inspect quality of reads.  This code only looks at two from forward and two from reverse as examples.
    #To view more, change the values.

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#e. Dada2 filter and trim process (I used default params) creates filtered files in filtered/ subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt_trimmed.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt_trimmed.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

#f. Dada2 learns error rates and plots them.  See tutorial for interpretation

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)


#g. Core sample dada2 algorithm inference

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#h. Merge reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#i. Inspect the merger data.frame from the first sample

head(mergers[[1]])

#j. Create sequence table of ASVs (seqtab), write as csv file.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
write.csv(seqtab, file = "seqtab.csv")

#k. Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

#l. Remove chimeras, and rewrite as new csv.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #For this run = .95874
write.csv(seqtab.nochim, file = "seqtab.nochim.csv")

#m. Track reads through dada2 pipeline, and create track.csv file (metadata) for phyloseq.  This 'track' file 
    #will need to have all your other metadata (from other sources) added to it before use downstream.

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, file = "track.csv")

#n. Assign Taxa to ASVs. I used the current Silva database.  Output is the taxa_table for phyloseq.

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa, file = "taxa.csv")

#o. Create an ASV table for future use (not required): Add header column 'ASV_No.' for later use in seqtab.nochim

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}


#p. Save an ASV fasta file

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

#q. Create and save ASV count tables

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA) #ASV table is a flipped version of the seqtab.nochim table
write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

##################################################################

#PART II: PHYLOSEQ - Create Phyloseq object from dada2 output tables (seqtab.nochim, taxa, and track), remove contaminants

library(Biostrings)
library(phyloseq)
library(2)
library(tidyverse)
library(dplyr)

packageVersion('phyloseq') #version 1.44 as of May 17 2023
theme_set(theme_bw())

#a. Read-in dada2 output tables if not already in global environment of R.

setwd("C:/Users/sj4879dj/OneDrive - MNSCU/MelchiorFiles/Research/Research - Sarrecenia Project/Sarr_2023/SarrMN_Final")

seqtab.nochim <- read.csv("seqtab.nochim.csv")
taxa <- read.csv("taxa.csv")
track <- read.csv("track.csv")

#b. Adjust file format 

row.names(seqtab.nochim) <- seqtab.nochim$ ASV
seqtab.nochim <- seqtab.nochim %>% select (- ASV)
seqtab.nochim <- as.matrix(seqtab.nochim)

row.names(taxa) <- taxa$ASV
taxa <- taxa %>% select (- ASV)
taxa <- as.matrix(taxa)

row.names(track) <- track$SampleID
track <- track %>% select (- SampleID)


#c. Name the three components for phyloseq creation, create phyloseq object; save as RDS file.

seqtab.nochim <- otu_table (seqtab.nochim, taxa_are_rows = T)
taxa <- tax_table(taxa)
track <- sample_data(track)

physeq <- phyloseq(seqtab.nochim, taxa, track)

saveRDS(physeq, file = "physeq.RDS")


##################################################################

#PART III.  Remove contaminant sequences with Decontam package

library(decontam)
packageVersion("decontam")

#a. Pull sample data into data frame format; plot library sizes for each sample

df <- as.data.frame(sample_data(physeq)) 
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize)) + geom_point()

#b. Identify contaminants by frequency (see documentation for algorithm)

library(decontam)

contamdf.freq <- isContaminant(physeq, method="frequency", conc="qPCR_log10")
head(contamdf.freq)
hist(contamdf.freq$p)
table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

#c. Plot Contaminants

plot_frequency(physeq, taxa_names(physeq)[c(1,4)], conc="qPCR_log10") + 
  xlab("DNA Concentration")

#d. Remove Contaminants

set.seed(100)

plot_frequency(physeq, taxa_names(physeq)[sample(which(contamdf.freq$contaminant),3)], conc="qPCR_log10") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

#e. New phyloseq object without contaminants

physeq.nc <- prune_taxa(!contamdf.freq$contaminant, physeq)
saveRDS(physeq.nc, "physeq_NC523.rds")



#f. Change ASV sequence to ASV1, ASV2, etc; create revised phyloseq object (complete)

library("vegan")
library(ggplot2)

dna <- Biostrings::DNAStringSet(taxa_names(PS1))
names(dna) <- taxa_names(PS1)
PS2 <- merge_phyloseq(PS1, dna)
taxa_names(PS2) <- paste0("ASV", seq(ntaxa(PS2)))
PS2
saveRDS(PS2, './physeq_NC523_ASVs.rds')


PS2=readRDS("physeq_NC523_ASVs.rds")
PS2

###################################################################################################

#PART IV:  Analysis of Alpha Diversity

#Read in ASV-level (total) phyloseq object
ps1 <- readRDS("physeq_NC_ASVs_alpha_LfYr703.rds")

##IV.1. Create phyloseq objects by taxon via agglomeration; add diversity indices, save as new phyloseq object.

#a. Agglomerate by taxon, calculate alpha diversity stats, add to metadata frames, and save new ps objects and rds files.

#Full (ASV-level) phyloseq object

#b (Calculate Pielou evenness with microbiome; add to ps object)

meta_ASV <- data.frame(sample_data(ps1))
alphaASV <-microbiome::alpha(ps1, index = "evenness_pielou")
head(alphaASV)


#c Calculate diversity matrices and add create final ASV-level ps object with all alpha div data.

meta_ASV <- data.frame(sample_data(ps_ASV.pie))
richASV <- estimate_richness(ps_ASV.pie, split = TRUE, measures = NULL)

PielouASV = alphaASV$evenness_pielou
meta_ASV$Pielou = PielouASV
Obs = richASV$Observed
meta_ASV$Obs = Obs
Chao1 = richASV$Chao1
meta_ASV$Chao1 = Chao1
Shannon = richASV$Shannon
meta_ASV$Shannon = Shannon
Simpson = richASV$Simpson
meta_ASV$Simpson = Simpson
InvSimpson=richASV$InvSimpson
meta_ASV$InvSimpson=InvSimpson
Fischer = richASV$Fisher
meta_ASV$Fischer = Fischer

meta_ASV = sample_data(meta_ASV)

physeq_NC_Rich626 = merge_phyloseq(ps_ASV.pie, meta_ASV)
physeq_NC_Rich626
saveRDS(physeq_NC_Rich626,"C:/Users/sj4879dj/OneDrive - MNSCU/MelchiorFiles/Research/Research - Sarrecenia Project/Sarr_2023/SarrMN_Final/physeq_NC_ASVs_alpha626.rds")


#d Calculate diversity matrices and add create final Family-level ps object with all alpha div data.

ps.fam <- tax_glom(ps1,"Family")

metafam <- data.frame(sample_data(ps.fam))
richfam <- estimate_richness(ps.fam, split = TRUE, measures = NULL)
piefam <-microbiome::alpha(ps.fam, index = "evenness_pielou")
Pielouf = piefam$evenness_pielou
metafam$Pielou.fam = Pielouf
Obs = richfam$Observed
metafam$Obs.fam = Obs
Chao1 = richfam$Chao1
metafam$Chao1.fam = Chao1
Shannon = richfam$Shannon
metafam$Shannon.fam = Shannon
Simpson = richfam$Simpson
metafam$Simpson.fam = Simpson
InvSimpson=richfam$InvSimpson
metafam$InvSimpso.fam=InvSimpson
Fischer = richfam$Fisher
metafam$Fischer.fam = Fischer

metafam= sample_data(metafam)
View(metafam)
ps.fam.1 = merge_phyloseq(ps.fam, metafam)
ps.fam.1

metacheck=data.frame(sample_data(ps.fam.1))
View(metacheck)

saveRDS(ps.fam.1, file = "C:/Users/sj4879dj/OneDrive - MNSCU/MelchiorFiles/Research/Research - Sarrecenia Project/Sarr_2023/SarrMN_Final/physeq.family_NC_alpha726.rds")


#IV.2  Testing alpha diversity metrics for distribution normality (QQ plots and Shapiro-Wilk test)

#a. Generate Q-Q Plots and Shapiro-Wilk Test for Normality

ggqqplot(meta$Obs) 
shapiro.test(meta$Obs) 

ggqqplot(meta$Pielou) 
shapiro.test(meta$Pielou) 

ggqqplot(meta$Shannon) 
shapiro.test(meta$Shannon) 


#IV.3  Mixed Effects Modeling of Alpha Diversity Metrics

      #Mixed-Effect Model Fitting.  Use commands from lme4 package to run mixed effects model. lmer model testing on all 
      #iterations of variables was accomplished and saved in LMER.code.fam.nlyr.R.

      #Chosen models:  The following were chosen because they included Sample_Date_No, Leaf Year (note: LeafYr, NOT LeafYr) AND Bog_No. as fixed effects, and LeafID as a 
      #Random effect. Names were changed from those above for more clarity.
      #Random effect (LeafID) generates variance = 0 (singular warnings), but does little to change fixed coefficients

library(lme4)
library(ggpubr)
library(lmerTest)
library(emmeans)

#a.  Read in ps object, convert sample data to data.frame (Note:This section uses an updated phyloseq object.  Name may differ from previous sections!!)

ps.fam.1=readRDS("physeq.family_NC_alpha726.rds")

metaASV=data.frame(sample_data(ps1))
metafam1=data.frame(sample_data(ps.fam.1))
View(metafam1)

#b. Final Models at ASV taxon level (model iteration testing code can be delivered upon request)

#Y =  Richness (Obs):  
Model.obs <- lmer(Obs ~ Sample_Date_No + + Bog_No. + LeafYr + qPCR_log10+ (1|LeafID), data=metaASV)
summary(Model.obs)
AIC(Model.obs)

#Y = Evenness (Pielou)
Model.pielou <- lmer(Pielou ~ Sample_Date_No + + Bog_No. + LeafYr + qPCR_log10+ (1|LeafID), data=metaASV)
summary(Model.pielou)
AIC(Model.pielou)

#Y = Shannon Diversity (Shannon):Model.shannon
Model.shannon <- lmer(Shannon ~ Sample_Date_No + + Bog_No. + LeafYr + qPCR_log10+ (1|LeafID), data=metaASV)
summary(Model.shannon)
AIC(Model.shannon)


#c.  Final Models at Family taxon level

#Y =  Richness (Obs):  
Model.fam.obs <- lmer(Obs.fam ~ Sample_Date_No + Bog_No. + LeafYr + qPCR_log10 + (1|LeafID), data=metafam1)
summary(Model.fam.obs)

#Y = Evenness (Pielou) (Sig vals for Intercept and Sept; AIC = -150)
Model.fam.pielou <- lmer(Pielou.fam ~ Sample_Date_No + Bog_No. + LeafYr + qPCR_log10 + (1|LeafID), data=metafam1) 
summary(Model.fam.pielou)
AIC(Model.fam.pielou)

#Y = Shannon Diversity (Shannon):Model.fam.shannon
Model.fam.shannon <- lmer(Shannon.fam ~ Sample_Date_No + Bog_No. + LeafYr + qPCR_log10 + (1|LeafID), data=metafam1)
summary(Model.fam.shannon)
AIC(Model.fam.shannon)

##In manuscript, Sample_Date_No, Bog_No., LeafYr, qPCR_log10, and LeafID are referred (for clarity and readability) 
##to Sample Date (Month), Bog, pitcher age, amplicon copies, and pitcher identity)

#AIC/BIC Testing
AIC(Model.fam.obs,Model.fam.pielou, Model.fam.shannon)
BIC(Model.fam.obs,Model.fam.pielou,Model.fam.shannon)

summary(Model.fam.obs)
summary(Model.fam.pielou)
summary(Model.fam.shannon)

#IV.4. Post Hoc Testing of Models for Significance

#a.  ANOVA (III) testing of models

anova3.Model.fam.obs <- anova(Model.fam.obs, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
anova3.Model.fam.obs

anova3.Model.fam.shannon <- anova(Model.fam.shannon, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
anova3.Model.fam.shannon

anova3.Model.fam.pielou <- anova(Model.fam.pielou, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
anova3.Model.fam.pielou


#b.  Eemeans and pair-wise testing

emmeans(Model.fam.obs, specs = pairwise~Sample_Date_No)
emmeans(Model.fam.obs, specs = pairwise~Bog_No.)

emmeans(Model.fam.shannon, specs = pairwise~Sample_Date_No)
emmeans(Model.fam.shannon, specs = pairwise~Bog_No.)

emmeans(Model.fam.pielou, specs = pairwise~Sample_Date_No)
emmeans(Model.fam.pielou, specs = pairwise~Bog_No.)

#################################
#Part V:  Plotting Alpha Diversity Metrics

library(ggplot2)
library(tidyverse)
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(scales)
library(ggpubr)

#a.  Generate histograms of alpha metrics to check distribution normality

Obs <- hist(metafam1$Obs.fam,ylim=c(0,20), breaks = 15)
Piel <- hist(metafam1$Pielou.fam ,ylim=c(0,30), breaks = 10)
Sh <- hist(metafam1$Shannon.fam, ylim=c(0,20), breaks = 10)

#b. Box Plot Graphs (Fig 2 a - f):  Alpha variable x sample date (Cumulative plot for all bogs, Family Level) with emmeans contrast p values

ps1 <- readRDS("physeq.family_NC_alpha726.rds")
ps1
ps1.fam=tax_glom(ps1,"Family")

metafam1 <- data.frame(sample_data(ps1.fam))

#b.1  Family richness by month

Obsxdate.fam <- ggplot(metafam1, aes(x=Sample_Date_No, y=Obs.fam)) + 
  geom_boxplot(outlier.shape = NA, color="blue4", fill="NA") + geom_jitter(position=position_jitter(0.2)) + xlab("Sample Month") +
  ylab("Richness (Families)") + 
  theme_classic() + scale_y_continuous(limits=c(0, 150))

Obsxdate.fam

Fig2a <- Obsxdate.fam + geom_bracket(xmin = "1 (May)", xmax = "4 (September)", y.position = 145,
                            label = "p < 0.01") + geom_bracket(xmin = "2 (July)", xmax = "4 (September)", y.position = 135,
                                                               label = "p < 0.001") + geom_bracket(xmin = "3 (August)", xmax = "4 (September)", y.position = 125,
                                                                                                label = "p < 0.001") + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0))

Fig2a

#Save as a high res eps vector graphic file:

ggsave("Fig2a2.eps", Fig2a, height = 2.5, width = 3, dpi = 600)

#b.2  Family richness by bog

Obsxbog.fam <- ggplot(metafam1, aes(x=Bog_No., y=Obs.fam)) + 
  geom_boxplot(outlier.shape = NA, colour = "red4", fill = "NA", ) + geom_jitter(position=position_jitter(0.2)) + xlab("Bog Location") +
  ylab("Richness (Families)")+theme_classic()+ scale_y_continuous(limits=c(0, 150))

Obsxbog.fam
Fig2d=Obsxbog.fam

ggsave("Fig2d.eps", Fig2d, height = 5, width = 6, dpi = 600)


#b.3  Shannon index (Family) by month

Shanxdate.fam <- ggplot(metafam1, aes(x=Sample_Date_No, y=Shannon.fam)) + 
  geom_boxplot(outlier.shape = NA, color="blue4", fill="NA") + geom_jitter(position=position_jitter(0.2)) + xlab("Sample Month") +
  ylab("Shannon Index") + ggtitle(label = "") + 
  theme_classic() + scale_y_continuous(limits=c(1, 5))

Shanxdate.fam

Fig2b <- Shanxdate.fam + geom_bracket(xmin = "1 (May)", xmax = "4 (September)", y.position = 4.6,
  label = "p < 0.001") + geom_bracket(xmin = "2 (July)", xmax = "4 (September)", y.position = 4.3,
  label = "p < 0.001") + geom_bracket(xmin = "3 (August)", xmax = "4 (September)", y.position = 4.0,
  label = "p < 0.001") +
  theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0))

Fig2b
ggsave("Fig2b.eps", Fig2b, height = 5, width = 6, dpi = 300)
ggsave("Fig2b.tiff", Fig2b, height = 5, width = 6, dpi = 600)

#b.4 Shannon index (family) by bog location

Shanxbog.fam <- ggplot(metafam1, aes(x=Bog_No., y=Shannon.fam)) + 
  geom_boxplot(outlier.shape = NA, colour = "red4", fill = "NA", ) + geom_jitter(position=position_jitter(0.2)) + xlab("Bog Location") +
  ylab("Shannon Index") +  theme_classic() + scale_y_continuous(limits=c(1, 5))
 
Fig2e=Shanxbog.fam

Fig2e
ggsave("Fig2e.eps", Fig2e, height = 5, width = 6, dpi = 300)

#b.5 Pielou evenness (Family) by month

Piexdate.fam <- ggplot(metafam1, aes(x=Sample_Date_No, y=Pielou.fam)) + 
  geom_boxplot(outlier.shape = NA, color="blue4", fill="NA") + geom_jitter(position=position_jitter(0.2)) + xlab("Sample Month") +
  ylab("Pielou Evenness Index") + theme_classic() + scale_y_continuous(limits=c(0.4, 0.9))

Piexdate.fam

Fig2c=Piexdate.fam + geom_bracket(xmin = "1 (May)", xmax = "4 (September)", y.position = 0.90,
                            label = "p < 0.01") + geom_bracket(xmin = "3 (August)", xmax = "4 (September)", y.position = 0.86,
                                                               label = "p < 0.001") +
  theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0))


Fig2c
ggsave("Fig2c.eps", Fig2c, height = 5, width = 6, dpi = 600)


#b.6 Pielou evenness (Family) by Bog

Piexbog.fam <- ggplot(metafam1, aes(x=Bog_No., y=Pielou.fam)) + 
  geom_boxplot(outlier.shape = NA, colour = "red4", fill = "NA", ) + geom_jitter(position=position_jitter(0.2)) + xlab("Bog Location") +
  ylab("Pielou Evenness Index") + 
  theme_classic() + scale_y_continuous(limits=c(0.4, 0.9))

Fig2f=Piexbog.fam

Fig2f
ggsave("Fig2f.eps", Fig2f, height = 5, width = 6, dpi = 600)


##PART V: Analysis of Beta Diversity

#Generate NMDS ordination plots (Bray-Curtis distances) by bog location (Bog_No.) and sample date/season (Sample_Date_No), then
#analyze dispersion of BC distances within/between bogs and dates (betadisper, PERMANOVA, ANOVA, and Tukey's post-hoc tests)

#V.1. Ordinate (NMDS) for Bray-Curtis at ASV and family level

library(phyloseq)
library(dplyr)

ps1=readRDS("physeq_NC_ASVs_alpha_LfYr703.rds")

psfam1=readRDS("physeq.family_NC_alpha726.rds")

#a.  Ordinate at ASV and Family levels

asv.ord <- ordinate(ps1, "NMDS", "bray", k=5)
fam.ord <- ordinate(psfam1, "NMDS", "bray", k=5)


#b. Plot NMDS Ordinations by bog number and sample date

library(tibble)
library(ggh4x)

#ASV level

asvord.date <- plot_ordination(ps1, asv.ord, type="samples", color="Sample_Date_No") + theme_classic() + theme(aspect.ratio = 1) + geom_point(size=2)+stat_ellipse() + scale_y_continuous(limits=c(-2.5, 2.2)) + scale_x_continuous(limits=c(-2.5, 2.2)) + force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))

asvord.date
Fig3a=asvord.date
ggsave("Fig3a.eps", Fig3a, height = 3, width = 5, dpi = 300)

asvord.bog <- plot_ordination(ps1, asv.ord, type="samples", color="Bog_No.") + theme_classic() + theme(aspect.ratio = 1) + geom_point(size=2) + stat_ellipse() + scale_y_continuous(limits=c(-2.5, 2.2)) + scale_x_continuous(limits=c(-2.5, 2.2)) + force_panelsizes(rows = unit(2, "in"), cols = unit(2, "in"))
Fig3b=asvord.bog
ggsave("Fig3b.eps", Fig3b, height = 3, width = 5, dpi = 300)


#Family level

famord.date <- plot_ordination(psfam704, fam.ord, type="samples", color="Sample_Date_No")+geom_point(size=3)+stat_ellipse()
famord.date

famord.bog <- plot_ordination(psfam704, fam.ord, type="samples", color="Bog_No.")+geom_point(size=3)+stat_ellipse()
famord.bog


#V.2 Statistical Analysis: Create Bray-Curtis distance matrices for each taxon (Bray-Curtis distance matrices) 

#a. Create matrices
library(vegan)

distasv = vegdist(t(otu_table(ps1)), method = "bray")
distasv

distfam = vegdist(t(otu_table(psfam1)), method ="bray" )
distfam

#b. Permutation ANOVA (PERMANOVA) via adonis2/vegan on modeled variables to partition (variation within the distance matrices?)

set.seed(1)

#Pull metadata (sample_data) from PS objects into data frame format:

meta_ASV <- data.frame(sample_data(ps1))
metafam <- data.frame(sample_data(psfam1))
View(meta_ASV)

#PERMANOVA

library(vegan)

permanova.asv_lf = adonis2(distasv ~ Sample_Date_No + Bog_No. + LeafYr + + PlantID + qPCR_log10 + LeafID, data = meta_ASV)
permanova.asv_lf

permanova.fam <- adonis2(distfam ~   Sample_Date_No + Bog_No. + LeafYr + PlantID + qPCR_log10 + LeafID, data = metafam)
permanova.fam


#c. Calculate mean dispersion around centroids by Sample Date and Bog No.

disp.month.asv = betadisper(distasv, type = c("centroid"), meta_ASV$Sample_Date_No)
disp.bog.asv = betadisper(distasv, type = c("centroid"), meta_ASV$Bog_No.)
disp.lf_yr.asv = betadisper(distasv, type = c("centroid"), meta_ASV$LeafYr)
disp.LeafID.asv = betadisper(distasv, type = c("centroid"), meta_ASV$LeafID)
disp.PlantID.asv = betadisper(distasv, type = c("centroid"), meta_ASV$PlantID)
disp.qPCR_log10.asv=betadisper(distasv, type = c("centroid"), meta_ASV$qPCR_log10)


disp.month.fam = betadisper(distfam, type = c("centroid"), metafam$Sample_Date_No)
disp.bog.fam = betadisper(distfam, type = c("centroid"), metafam$Bog_No.)
disp.lf_yr.fam = betadisper(distasv, type = c("centroid"), metafam$LeafYr)
disp.LeafID.fam = betadisper(distasv, type = c("centroid"), metafam$LeafID)
disp.PlantID.fam = betadisper(distasv, type = c("centroid"), metafam$PlantID)
disp.qPCR_log10.fam=betadisper(distasv, type = c("centroid"), metafam$qPCR_log10)

#d. ANOVAs on dispersion by Sample Date and Bog No.

anova(disp.month.asv)
anova(disp.bog.asv)
anova(disp.lf_yr.asv)
anova(disp.LeafID.asv)
anova(disp.PlantID.asv)
anova(disp.qPCR_log10.asv)

anova(disp.month.fam)
anova(disp.bog.fam)
anova(disp.lf_yr.fam)
anova(disp.LeafID.fam)
anova(disp.PlantID.fam)
anova(disp.qPCR_log10.fam)


#e. Tukey's test (post-hoc) to determine where sig difs are (if found by ANOVA) by Sample Date and Bog

TukeyHSD(disp.LeafID.asv, which = "group", ordered = FALSE,
         conf.level = 0.95)

TukeyHSD(disp.PlantID.asv, which = "group", ordered = FALSE,
         conf.level = 0.95)

TukeyHSD(disp.month.asv, which = "group", ordered = FALSE,
         conf.level = 0.95)

TukeyHSD(disp.bog.asv, which = "group", ordered = FALSE,
         conf.level = 0.95)

TukeyHSD(disp.month.fam, which = "group", ordered = FALSE,
         conf.level = 0.95)

TukeyHSD(disp.bog.fam, which = "group", ordered = FALSE,
         conf.level = 0.95)

TukeyHSD(disp.lf_yr.asv, which = "group", ordered = FALSE,
         conf.level = 0.95)


###################################################################################################

##PART VI:  Differential Abundance Testing and Plots

#VI.1  CSS normalization (metagenomeSeq package) of count data for diff abundance testing; save new phyloseq object

library(vegan)
library(metagenomeSeq)
library(Maaslin2)
library(dplyr)

ps1 <- readRDS("physeq_NC_ASVs_alpha_LfYr703.rds")

#a. Normalize data and create new phyloseq object 

metaseq_asv_tab <- phyloseq_to_metagenomeSeq(ps1)
metaseq_asv_tab_css <- cumNorm(metaseq_asv_tab, p=cumNormStatFast(metaseq_asv_tab))
metaseq_asv_tab_css=MRcounts(metaseq_asv_tab_css, norm=T)

PS_CSS2 <- merge_phyloseq(otu_table(metaseq_asv_tab_css,
                                    taxa_are_rows=T), tax_table(ps1), sample_data(ps1))

saveRDS(PS_CSS2, './physeq_NC703_ASV_CSS.rds')
PS_CSS2  ## Use this ps object only for diff abundance testing.

PS_CSS2= readRDS('./physeq_NC703_ASV_CSS.rds')

#b.  Agglomerate normalized count data at the family level, and create a new ps object

PS_CSS2_fam=tax_glom(PS_CSS2, "Family")

#VI.2.  Differential Abundance Testing with Maaslin2 (family level)

#a.  Reformatting data for Maaslin2

#Create input asv file from ps otu table:
input_data_asv=data.frame(otu_table(PS_CSS2_fam))

#Create input metadata file from ps sample data table:
input_metadata = data.frame(sample_data(PS_CSS2_fam))

#Set Bog (primary fixed effect variable) as a factor 
input_metadata$Bog = factor(input_metadata$Bog, levels=c("BBSR_North", "Clearwater_Mid", "Beckman_South"))

#Create input tax table file from ps object:
tax_table_fam=tax_table(PS_CSS2_fam)[,"Family"] 

#Create and add new variable column 'ASVID' to the tax table, and set table as a data frame.
tax_table_fam= data.frame(tax_table_fam)%>%rownames_to_column(var='ASVID')

input_data_asv_edit=input_data_asv%>% 
  rownames_to_column(var='ASVID')

#Adding (joining) ASV name data to tax table
input_data_asv_edit_left_join=left_join(input_data_asv_edit, tax_table_fam, by=c("ASVID"))

input_data_asv_edit_left_join$ASVID <-input_data_asv_edit_left_join$Family
colnames(input_data_asv_edit_left_join)
df = subset(input_data_asv_edit_left_join, select = -c(Family))
input_data_asv_edit_df=df%>%column_to_rownames(var='ASVID')
input_data_asv_edit_df_tranpose=t(input_data_asv_edit_df)
 
#c.  Run Maaslin2: Input files (Family count table and metadata). 

#Fixed effect under consideration is Bog

fit_data2 = Maaslin2(
  input_data = input_data_asv_edit_df_tranpose, 
  input_metadata = input_metadata, 
  output = "Bog_BBSR_Ref", 
  fixed_effects = c("Bog"),
  random_effects = c("LeafID"),
  min_prevalence = 0.1,
  min_abundance = 0.01)

#Fixed effect under consideration is Sample_Date_No

#Set Sample_Date_No (primary fixed effect variable) as a factor 
input_metadata$Sample_Date_No = factor(input_metadata$Sample_Date_No, levels=c("1 (May)","2 (July)","3 (August)", "4 (September)"))

fit_data2 = Maaslin2(
  input_data = input_data_asv_edit_df_tranpose,
  input_metadata = input_metadata,
  output = "Date_May_Ref",
  fixed_effects = c("Sample_Date_No"),
  random_effects = c("LeafID"),
  min_prevalence = 0.1,
  min_abundance = 0.01)



#VI.3.  Differential Abundance Testing with Maaslin2 (genus level)

#a.  Agglomerate normalized count data at the genus level, and create a new ps object

css_gen.ps=tax_glom(PS_CSS2, "Genus")
variable.names(tax_table(css_gen.ps))
css_gen.is.ps <- subset_taxa(css_gen.ps, Genus != "Incertae Sedis")

#b.  Reformatting data for Maaslin2 (genus level)

#Create input ASV file from ps OTU table:
input_data_asv=data.frame(otu_table(css_gen.is.ps))

#Create input metadata file from ps sample data table:
input_metadata = data.frame(sample_data(css_gen.is.ps))

#Set Sample_Date_No (primary fixed effect variable) as a factor 
input_metadata$Sample_Date_No=as.factor(input_metadata$Sample_Date_No)
input_metadata$Sample_Date_No = factor(input_metadata$Sample_Date_No, levels=c("1 (May)","2 (July)","3 (August)", "4 (September)"))

#Create input tax table file from ps object:
tax_table_gen=tax_table(css_gen.is.ps)[,"Genus"] 

#Create and add new variable column 'ASVID' to the tax table, and set table as a data frame.
tax_table_gen= data.frame(tax_table_gen) %>% rownames_to_column(var='ASVID')
input_data_asv_edit=input_data_asv %>% rownames_to_column(var='ASVID')


#Adding (joining) ASV name data to tax table
input_data_asv_edit_left_join=left_join(input_data_asv_edit, tax_table_gen, by=c("ASVID"))

input_data_asv_edit_left_join$ASVID <-input_data_asv_edit_left_join$Genus
colnames(input_data_asv_edit_left_join)

df = subset(input_data_asv_edit_left_join, select = -c(Genus))
input_data_asv_edit_df = df %>% column_to_rownames(var='ASVID')
input_data_asv_edit_df_tranpose=t(input_data_asv_edit_df)

#d. Run Maaslin2:  Input files (Genus count table and metadata). 

#Fixed effect under consideration is Sample_Date_No

fit_data2 = Maaslin2(
  input_data = input_data_asv_edit_df_tranpose, 
  input_metadata = input_metadata, 
  output = "Genus_month", 
  fixed_effects = c("Sample_Date_No"),
  random_effects = c("LeafID"),
  min_prevalence = 0.1,
  min_abundance = 0.01)


#Fixed effect under consideration is Sample_Date_No

#Create input asv file from ps otu table:
input_data_asv=data.frame(otu_table(css_gen.is.ps))

#Create input metadata file from ps sample data table:
input_metadata = data.frame(sample_data(css_gen.is.ps))

#Set Sample_Date_No (primary fixed effect variable) as a factor 
input_metadata$Bog=as.factor(input_metadata$Bog)
input_metadata$Bog = factor(input_metadata$Bog, levels=c("Beckman_South", "Clearwater_Mid", "BBSR_North"))

#Create input tax table file from ps object:
tax_table_gen=tax_table(css_gen.is.ps)[,"Genus"] 

#Create and add new variable column 'ASVID' to the tax table, and set table as a data frame.
tax_table_gen= data.frame(tax_table_gen) %>% rownames_to_column(var='ASVID')

input_data_asv_edit=input_data_asv %>% rownames_to_column(var='ASVID')

input_data_asv_edit_left_join=left_join(input_data_asv_edit, tax_table_gen, by=c("ASVID"))

input_data_asv_edit_left_join$ASVID <-input_data_asv_edit_left_join$Genus
colnames(input_data_asv_edit_left_join)

df = subset(input_data_asv_edit_left_join, select = -c(Genus))

input_data_asv_edit_df = df %>% column_to_rownames(var='ASVID')

input_data_asv_edit_df_tranpose=t(input_data_asv_edit_df)


fit_data2 = Maaslin2(
  input_data = input_data_asv_edit_df_tranpose, 
  input_metadata = input_metadata, 
  output = "Genus_Bog", 
  fixed_effects = c("Bog"),
  random_effects = c("LeafID"),
  min_prevalence = 0.1,
  min_abundance = 0.01)


#VI.4:  Plotting Differential Abundance Results

library(lemon)          
library(ggrepel)
library(ggplot2)

#a.  Read in data from edited Maaslin2 data with RA means and prev added:

dmo <- read.csv("DiffAbund_by_Month_714.csv")

#b. Create header and da objects (title and da for plotting)

header = "(a) Differential Abundance By Month (Pairwise Comparisons to May)"
da = "May"

#d. Create new columns for 'Significant' q values, and for rel abundance in various ranges

dmo$Significant <- ifelse(dmo$q.value < 0.05, "q.value < 0.05", "Not Significant")
View(dmo)
dmo$Mean01.RA <- ifelse(dmo$Mean.RA < .1, "RA < .1", "RA > .1")  ## Here, the RAs in data are as percentages, so '0.1' = 0.1% RA

#dmo$Mean05.RA <- ifelse(dmo$Mean.RA < .5, "RA < .5", "RA > .5")
#dmo$Mean1.RA <- ifelse(dmo$Mean.RA < 1, "RA < 1", "RA > 1")


#e. Filter in data over RA thresholds.  These will be plotted.  

dmo01 <- subset(dmo, dmo$Mean.RA > .1)

#f. Set facet grouping

dmo$group <- factor(dmo$group,levels=c("July","August","September"))

#g.  Differential relative abundance (Family) plot (Fig 6a and 6b) by month 

bplot01 <- ggplot(dmo01, aes(x = Log2, y = Family, fill = Significant)) +
  #geom_point(aes(size=size_AveEx),width = 0.30, height = 0.35, alpha = 1, na.rm = T, shape = 21, colour = "black") +
  geom_point(aes(size=Mean.RA), alpha = 1, na.rm = T, shape = 21, colour = "black") +
  scale_fill_manual(values = c("grey70","red"))+       
  theme(
    panel.grid.major = element_line(colour="#f0f0f0"),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=10, margin = margin(.5,0,.5,0)),
    strip.text.y=element_text(size=10, margin = margin(.5,0,.5,0)),
    axis.text.y=element_text(size=10),
    axis.title=element_text(size=10),
    legend.position="right",
    panel.spacing=unit(0.3, "lines"),
    plot.title=element_text(size=12, hjust=0.5),
    legend.text=element_text(size=8),
    legend.title=element_text(size=10),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA,color = "grey")) + 
  xlim(-4,4)+
  geom_vline(xintercept = 0, linetype="dashed")+
  #facet_wrap(~Month, nrow=1)+  #Used if the 'Months' in the table are in the desired order
  facet_wrap(~factor(Month, levels=c('July','August','September'))) +  #Reorder facets
  ylab("Family")+xlab("Log2 Abundance Differential")+
  theme(legend.position = 'none')


bplot01
ggsave("Fig6a.eps", bplot01, height = 6, width = 6, dpi = 600)

#h.  Differential relative abundance (Family) plots (Fig 6a and 6b) by month 

#i. Read in data from edited Maaslin2 data with RA means and prev added:
dbog2 = read.csv("DiffAbund_bog_714_Adj.csv")

#j.  Create a title (header) and da objects for plotting
headerbog = ""

#k. Create new columns for 'Significant' q values, and for rel abundance in various ranges, and 
#data columns for different mean rel abundance cuttoffs

dbog2$Significant <- ifelse(dbog2$q.value < 0.05, "q.value < 0.05", "Not Significant")
dbog2$Mean01.RA <- ifelse(dbog$Mean.RA < .1, "RA < .01", "RA > .01") 

#l. Filter dataframe for taxa with RA > 'X' values.  These will be plotted by ggplot2

dbog201 <- subset(dbog2, dbog2$Mean.RA > .1

#m.  Differential relative abundance (Family) plot (Fig 6a and 6b) by month                 
                  
dbplot201 <- ggplot(dbog201, aes(x = Log2.Differential, y = Family, fill = Significant))+ 
  geom_point(aes(size=Mean.RA), alpha = 1, na.rm = T, shape = 21, colour = "black") +
  scale_fill_manual(values = c("grey70", "red"))+
  theme(
    panel.grid.major = element_line(colour="#f0f0f0"),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=10, margin = margin(.5,0,.5,0)),
    strip.text.y=element_text(size=10, margin = margin(.5,0,.5,0)),
    axis.text.y=element_text(size=10),
    axis.title=element_text(size=10),
    legend.position="right",
    panel.spacing=unit(0.3, "lines"),
    plot.title=element_text(size=12, hjust=0.5),
    legend.text=element_text(size=8),
    legend.title=element_text(size=10),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA,color = "grey")) + 
  xlim(-4,4)+
  geom_vline(xintercept = 0, linetype="dashed")+
  facet_wrap(~Ref_Bog.Comp_Bog, nrow=1)+  #This creates multi panels of plots
  ylab("Family")+xlab("Log2 Abundance Differential")+
  ggtitle(headerbog)+         
  theme(legend.position = 'none')

dbplot201
ggsave("Fig6b.eps", dbplot201, height = 6, width = 6, dpi = 600)


#n. Volcano Plot (by Month) - Genera

dgmo <- read.csv("DA.genus.month719.csv")

#o. Create header and da objects (title and da for plotting)

header2 = "Differential Genus Abundance By Month (Pairwise Comparisons to May)"


#p. Create new columns for 'Significant' q values, and for rel abundance in various ranges
View(dgmo)
dgmo$Significant <- ifelse(dgmo$q.value < 0.05, "q.value < 0.05", "Not Significant")
dgmo$Mean01.RA <- ifelse(dgmo$Mean.RA < .01, "RA < .01", "RA > .01")

#q. Filter in data over RA thresholds, set values.  These will be plotted.  

dgmo01 <- subset(dgmo, dgmo$Mean.RA > .01)
View(dgmo01)

#Set relative abundance values as numeric
dgmo01$Rel.Abund<-as.numeric(dgmo01$Rel.Abund)
dgmo01$Mean.RA<-as.numeric(dgmo01$Mean.RA)

#Change the NA values to 0
dgmo01$Mean.RA[is.na(dgmo01$Mean.RA)] <- 0
dgmo01$Rel.Abund[is.na(dgmo01$Rel.Abund)] <- 0

#Check structure of variables:
str(dgmo01$Mean.RA) #Etc.

# Re-leveling categorical variable for plotting (order of plots)
relevel=c("July","August","September")
dgmo01$Month = factor(dgmo01$Month, levels=relevel)

#r.  Plot

DA.Gen.volcano <- ggplot(dgmo01, aes(x=Log2.Differential, y=-log10(q.value), fill=Significant))+
  geom_point(aes(size=Mean.RA), alpha=1, shape=21)+
  scale_fill_manual(values=c( "grey87","red"))+
  geom_vline(xintercept = 0, linetype="dashed")+
  ggtitle(header2)+
  labs(fill="Significance", size="Mean Rel. Abund. (%)")+
  #ylim(0.5,4.0)+
  facet_wrap(~Month)+ # you can activate this if you want graphs facetted
  geom_text_repel(data=subset(dgmo01, q.value < 0.05), aes(label=Genus), size=2,box.padding=unit(0.2, "lines"), point.padding=unit(0.2, "lines"))

DA.Gen.volcano


###############################################################################################################

#PART VIII:  Relative Abundance and Prevalence Plots

library(tidyverse)
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(scales)
library(paletteer)

#VIII.1:  Set up data for plotting

#a.  Read in and agglomerate non-normalized phyloseq data at family and genus levels.  Create new PS objects

psfull=readRDS("physeq_NC_ASVs_alpha_LfYr703.rds")
psfam=tax_glom(psfull, "Family", NArm = F)
psgen=tax_glom(psfull, "Genus", NArm = F)

#b.  Create plotting groups (by Bog and by Sample Date)
sample_data(psfam)$Sample_Date_No <- factor(sample_data(psfam)$Sample_Date_No)
sample_data(psfam)$Bog_No. <- factor(sample_data(psfam)$Bog_No.)

sample_data(psgen)$Sample_Date_No <- factor(sample_data(psgen)$Sample_Date_No)
sample_data(psgen)$Bog_No. <- factor(sample_data(psgen)$Bog_No.)

#VIII.2: Relative Abundance (Compositional) by Family 

#a. Generate mean RA by family

# Start with original ASV ps object; agglom to family, transform counts to RA x 100 (to get %)
# The taxa orders are the same in tax_table() and otu_table(); 'NArm = F' to keep NA counts in.
# Use 'rowMeans(otu_table(psobject))' to calculate per-row average in the OTU table

famform = psfam %>% transform_sample_counts(function(x) {x * 100/sum(x)})
dff = data.frame(Family = tax_table(famform)[,"Family"], Mean = rowMeans(otu_table(famform)), row.names = NULL)
dff = dff[order(-dff$Mean),]
write.csv(dff,"MeanRAs_fam_1223.csv")

#b. Generate mean RA by genus

genform = psgen %>% transform_sample_counts(function(x) {x * 100/sum(x)})
dfg = data.frame(Genus = tax_table(genform)[,"Genus"], Mean = rowMeans(otu_table(genform)), row.names = NULL)
dfg = dfg[order(-dfg$Mean),]
dfg=write.csv(dfg, "MeanRAs_genus_1223.csv")

read.csv("MeanRAs_fam_1223.csv")
read.csv("MeanRAs_genus_1223.csv")

#c.  Use microbiome pkg to transform to RA for both family and genus levels

ps.fam.RA <- psfam %>% aggregate_taxa(level = "Family") %>%
  microbiome::transform(transform = "compositional")
ps.fam.RA

ps.genus.RA <- psgen %>% aggregate_taxa(level = "Genus") %>%
  microbiome::transform(transform = "compositional")
ps.genus.RA

#d. Aggregate top (e.g., n = 20) families by relative abundance (note:  mean RAs are the same when calculated by the 
#phyloseq command and the microbiome commands.  Ref is 'microbacteriaceae at 9.8375 for both)

ps.fam.RA20 <- ps.fam.RA %>% aggregate_top_taxa2(level = "Family", top=20) %>%
  microbiome::transform(transform = "compositional")

fammic=as.data.frame(otu_table(ps.fam.RA))
fammic=write.csv(fammic, "Family_RAs_717.csv")

#Same, but for top 200 to get all for table analysis

ps.fam.RA200 <- ps.fam.RA %>% aggregate_top_taxa2(level = "Family", top=200) %>%
  microbiome::transform(transform = "compositional")

fammic2=as.data.frame(otu_table(ps.fam.RA200))
fammic2=write.csv(fammic, "Top200fams_RA.csv")

#VIII.3: Plotting Relative abundance 

#a.  Set up color palette for plot; use Paletteer package in R

Paul1 <- paletteer_d("colorBlindness::SteppedSequential5Steps") 

#b. Plot Relative Abundance (y) by sample date (x1) and/or bog location (x2)

Fig4a=ps.fam.RA20 %>%  plot_composition(group_by = "Sample_Date_No")+ scale_y_continuous(labels = percent)+ scale_fill_manual(values = Paul1)
Fig4a
ggsave("Fig4a.eps", Fig4a, height = 4, width = 9.5, dpi = 600)

Fig4b=ps.fam.RA20 %>%  plot_composition(group_by = "Bog_No.")+ scale_y_continuous(labels = percent)+ scale_fill_manual(values = Paul1)
Fig4b
ggsave("Fig4b.eps", Fig4b, height = 4, width = 9.5, dpi = 600)


#c.  Plot averages by bog and sample date

ps.fam.RA20 %>%  plot_composition(average_by = "Bog_No.")+ scale_y_continuous(labels = percent)
ps.fam.RA20 %>%  plot_composition(average_by = "Sample_Date_No")+ scale_y_continuous(labels = percent)


#VIII.4:  Calculate and plot taxon prevalence

#a. Calculate prevalence of families in samples

prev.fam30 <- prevalence(ps.fam.RA30,detection=0,sort=TRUE,count=FALSE)
prev.fam30=write.csv(prev.fam30, file="prevfam30.csv")
frame30 <- data.frame(read.csv("prevfam30.csv"))

#b. Plot prevalence (top 30, across all samples) - families

plotfam30 <- ggplot(data=frame30, aes(x=reorder(Family, + Prevalence), y=Prevalence, fill=Prevalence)) +
  geom_bar(stat="identity", fill="chocolate3") +
  coord_flip() +
  ylim(0, 1) +
  theme_classic()

plotfam30

plotfam30l=plotfam30 + labs(y = "Prevalence", x = "Family")
plotfam30l

#c.  Plot prevalence - diazotrophic genera

#(Di.gent is sheet with >2% prev listed only)

Di.gen <- data.frame(read.csv("Diazotrophs_RA_Genera_Sarr.csv"))
Di.gent <- data.frame(read.csv("Diazotrophs_RA_Genera_Sarr_T30.csv"))

Prev.diaz.gen <- ggplot(data=Di.gent, aes(x=reorder(Genus, +Prevalence), y=Prevalence, fill = Prevalence)) +
  geom_bar(stat="identity", fill='darkgreen')+
  coord_flip()+
  ylim(0, 1) +
  theme_classic()
Prev.diaz.gen
Prev.diaz.gent=Prev.diaz.gen + labs(y = "Prevalence", x = "Diazotrophic Genera")
Prev.diaz.gent

#Save as a high resolution eps file

Fig5=Prev.diaz.gent
ggsave("Fig5.eps", Fig5, height = 4.6, width = 7, dpi = 300)


################################################################################################################
#PART VII:  Analysis of qPCR data for nifH gene quantification

library(lme4)
library(emmeans)

##Used Excel to add in adjusted/normalized nifH quantification values, and exported as 'meta9.csv' to R

#a. Model relationships between date, bog, etc vs.nifH

meta10 <- read.csv("meta9_final.csv")

Model.nif <- lmer(Log10.nifH ~ Sample_Date_No + LeafYr + qPCR_log10 + Bog_No. + (1|LeafID), data=meta10)
summary(Model.nif)
AIC(Model.nif)  #AIC for this model is 253.115

#b.  Type III ANOVA testing of model

anova3.Model.nif <- anova(Model.nif, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
anova3.Model.nif

#c.  eemeans and pair-wise testing of contrasts

emmeans(Model.nif, specs = pairwise~Sample_Date_No)
emmeans(Model.nif, specs = pairwise~Bog_No.)


#d. Boxplot of nifH x sample date

meta10=read.csv("meta9_final.csv")
as.data.frame(meta10)


##nifH x sample date

nif.plot.mo <- ggplot(meta10, aes(x=Sample_Date_No, y=Log10.nifH)) + geom_boxplot(outlier.shape = NA, color="darkgreen", fill="NA") + geom_jitter(position=position_jitter(0.2)) + xlab("Sample Month") +
  ylab("Log10 nifH copies/mL") + 
  theme_classic()

nif.plot.mo
Fig7a=nif.plot.mo + geom_bracket(xmin = "1 (May)", xmax = "3 (August)", y.position = 7.75,
                                  label = "p < 0.001") + geom_bracket(xmin = "2 (July)", xmax = "3 (August)", y.position = 7.40,
                                                                     label = "p < 0.001") + geom_bracket(xmin = "3 (August)", xmax = "4 (September)", y.position = 7.05,
                                                                                                         label = "p < 0.01") +
  theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0))


ggsave("Fig7a.eps", Fig7a, height = 5, width = 6, dpi = 600)


##nifH x bog

nif.plot.bog <- ggplot(meta10, aes(x=Bog_No., y=Log10.nifH)) + 
  geom_boxplot(outlier.shape = NA, color="darkgreen", fill="NA") + geom_jitter(position=position_jitter(0.2)) + xlab("Sample Month") +
  ylab("Log10 nifH copies/mL") + 
  theme_classic()

nif.plot.bog
ggsave("Fig7b.eps", nif.plot.bog, height = 5, width = 6, dpi = 300)




