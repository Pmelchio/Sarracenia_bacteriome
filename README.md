Analysis of the Northern Pitcher Plant (Sarracenia purpurea L.) phytotelm bacteriome throughout a temperate region growing season 
Paul Melchior (C) 2024

In this investigation, we used high-throughput 16S rDNA sequencing and bioinformatics to analyze the S. purpurea phytotelm bacteriome at different time points through 
the growing season (May – September) in plants from the north-central region of the species’ native range (Minnesota). Additionally, qPCR to detect and quantify bacterial 
nitrogenase genes (nifH) in all phytotelm samples.  All bioinformatic and statistical analysis was done in R.

Table of R Code Contents (R packages used)

PART I:  	16S rDNA sequence filtering, trimming, and pair-end read joining, and taxonomy assignments

  		I.1.  Removal of primers (cutadapt)
  		I.2.  16S sequence filteration, additional trimming, PE read joining, chimera removal, and taxonomy (dada2)

PART II: 	Create metadata, sequence and taxa tables (phyloseq)

PART III:  	Revove contaminant sequences (decontam)

PART IV:	Alpha diversity analysis 

		IV.1.  Create phyloseq objects, calculate diversity matrices and add to PS objects (phyloseq)
		IV.2.  Test alpha diversity metrics for distribution normality (phyloseq)
		IV.3.  Mixed Effect Modeling of alpha metrics (lme4, lmerTest, emmeans)
		IV.4.  Post-hoc Testing (ANOVA 3; Tukey, est. marginal means)

PART V:		Plotting alpha diversity - Manuscript Fig 2 (ggplot2, microbiome, microbiomeutilities, scales, and ggpubr)

PART VI:	Beta diversity analysis - Manuscript Fig 3

		VII.1. Ordinate (NMDS) for Bray-Curtis Inices (phyloseq)
		VII.2. Creat Bray-Curtis distance matrices, PERMANOVA, and beta dispersion and post-hoc analysis (vegan)

PART VII.	Differential abundance testing and plots - Manuscript Fig 6

		VII.1. CSS Normalization (metagenomeSeq, phyloseq)
		VII.2. Differential abundance testing (maaslin2) - family level
		VII.3. Differential abundance testing (maaslin2) - genus level
		VII.4. Plotting differential abundance (ggplot2, ggrepel, lemon)

PART VIII:	Relative abundance and prevalence - Manuscript Fig 4 and Fig 5

		VIII.1. Setting up data for plotting (phyloseq)
		VIII.2. Generate RA by family (phyloseq, microbiome, microbiomeutilities)
		VIII.3. Plotting RA (ggplot2, paletteer, scales)
		VIII.4. Calculate and plot taxon prevalence

PART IX:	Analysis of qPCR data for nifH gene quantification - Manuscript Fig 7 (lme4, emmeans)
