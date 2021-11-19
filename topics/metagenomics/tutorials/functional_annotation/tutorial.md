---
layout: tutorial_hands_on

title: Preprocessing for microbiota shotgun data
time_estimation: 60m
questions:
- 
objectives:
- 
key_points:
- 
requirements:
contributors:
- bebatut
---

*What are They Capable of Doing?* Inferring Biological Function

![](images/functional_annotation.png "Figure 3 from [Sharpton, 2014](https://www.frontiersin.org/articles/10.3389/fpls.2014.00209/full) & Figure 1 from [Kim et al, 2013](https://www.researchgate.net/publication/257757057_Kim_M_Lee_KH_Yoon_SW_Kim_BS_Chun_J_Yi_H_Analytical_tools_and_databases_for_metagenomics_in_the_next-generation_sequencing_era_Genomics_Inform_11_102-113)")


genetic repertoire of a microbial community can be identified using adapted single-genome characterization tools

1. gene identification step, usually with a metagenomic-specific parameter setting
2. homology-based annotation pipelines commonly used for characterizing pure isolate genome assemblies

[Niu et al 2017](https://academic.oup.com/bib/article/3805128)


# Gene prediction

1. Gene fragment recruitment ~ Supervised binning
3. *de novo* gene prediction
    - determines which metagenomic reads contain coding sequences
    - Statistical methods using GC content, di-coding frequencies or translation table, hexamers statistics, ribosomal binding side motifs
        - Logistic regression models 
            - MetaGene ([Noguchi et al, 2006](https://academic.oup.com/nar/article-abstract/34/19/5623/3112023))
            - MetaGeneAnnotator ([Noguchi et al, 2008](https://academic.oup.com/dnaresearch/article/15/6/387/512877))
        - Log-likelihood function
            - MetaProdigal ([Hyatt et al, 2012](https://academic.oup.com/bioinformatics/article/28/17/2223/246063))

                Using of a log-likelihood function to distinguish the coding regions from the background

    - Model-based methods
        - HMM
            - MetaGeneMark ([Zhu ey al, 2010](https://academic.oup.com/nar/article/38/12/e132/2409881))
                1. Training genomes used to estimate a polynomial and logistic approximations of oligonucleotides freq as function of GC-content
                2. Estimates used as parameters of HMM
            - FragGeneSCan ([Rho et al, 2010](https://academic.oup.com/nar/article/38/20/e191/1317565))

                integration of codon usage bias, start/stop codon patterns, explicit modeling of sequencing errors

        - IMM
            - Glimmer-MG ([Kelley et al, 2012](https://academic.oup.com/nar/article/40/1/e9/1270026))
        - Neural network
            - Orphelia ([Hoff et al, 2009](https://academic.oup.com/nar/article/37/suppl_2/W101/1134091))
            

# Functional annotation

> **Protein family**: group of sequences that share a common evolutionary origin, reflected by their related functions and similarities in sequence or structure
> 
> **Gene ontology**: a major bioinformatics initiative to unify the representation of gene and gene product attributes across all species (controlled vocabulary, etc)
> 
> **Pathway**: a series of actions among molecules in a cell that leads to a certain product or a change in a cell

![](images/functional_annotation_db.png "Figure 1 from [Kim et al, 2013](https://www.researchgate.net/publication/257757057_Kim_M_Lee_KH_Yoon_SW_Kim_BS_Chun_J_Yi_H_Analytical_tools_and_databases_for_metagenomics_in_the_next-generation_sequencing_era_Genomics_Inform_11_102-113)")


Characterization of protein family: comparison of full-length protein sequences identified through genome sequencing projects

Similar biological functions

IMG/MER utilizes HMMsearch (profile HMMs) to associate genes with PFAM, and genes are further annotated using COGs. Database of position-specific scoring matrix (PSSMs) for COGs are downloaded from NCBI and are used to annotate protein sequences. Moreover, genes are labeled using KEGG-associated KO terms, EC numbers, and assigned phylogeny using similarity searches. With a large set of genomes in its public repositories, IMG/MER can exploit its own resources, using them as reference nonredundant databases from which it obtains additional functional annotation.



![](images/humann2.jpg "[HUMAnN2](https://bitbucket.org/biobakery/humann2/wiki/Home)")


![](images/humann2_pathway.png "Figure 3 from [Abubucker et al, 2012](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002358)")

Metabolic modules differentially present or abundant in at least one body habitat of the human microbiome