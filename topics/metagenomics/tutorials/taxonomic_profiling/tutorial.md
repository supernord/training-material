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


# Introduction

**Taxonomic profiling**: estimation of the abundances of the different taxa found in a metagenomics sample (not the classification of each read)

![](images/taxonomic_profile_example.png)

Applications
- Quantitative community profiling
- Clinical microbiology ([Schlaberg et al, 2017](http://www.archivesofpathology.org/doi/abs/10.5858/arpa.2016-0539-RA?code=coap-site), [Salzberg et al, 2016](http://nn.neurology.org/content/3/4/e251.short))
    - Presence or absence of infectious pathogens, which can be identified by matching reads against a reference database
    - Even though human-associated microbes are comparatively well studied with many completed genomes in the reference database, some pathogens remain unsequenced, and others have only recently been discovered using metagenomics sequencing


# Marker gene analysis

**Orthologous genes**: genes in different species that originated by vertical descent from a single gene of the last common ancestor

![](images/Ortholog_paralog_analog_examples.svg "Image by [Thomas Shafee](https://commons.wikimedia.org/w/index.php?curid=68505353)")

**Marker gene**: an orthologous gene group which can be used to delineate between taxonomic lineages

![](images/marker_gene_analyis.png "Figure 2 from [Sharpton, 2014](https://www.frontiersin.org/articles/10.3389/fpls.2014.00209/full)")

???
Marker gene approaches: identification of sets of clade-specific, single-copy genes (identification of one of these genes can be used as evidence that a member of the associated clade is present)

If a marker happens to be part of a fragment --> the entire fragment can be anchored by that marker to a specific taxonomic clade

*Which genes is a good taxonomic marker?*

???
16S rRNA

---

### Marker gene analysis

![](images/PhyEco_markers.png) "Figure 3 from [Wu et al, 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0077033)")

Which genes should be there and not be there to?

- Identify an Archaea
- Identify a Epsilonproteobacteria
- Distinguish a Beta-gamma-proteobacteria from an Epsilonproteobacteria


- Identify Archaea: none of the markers
- Identify Epsilonproteobacteria: everything should be there except B000064, B000076, B000087, B000092, B000098, B000106, B000112, B00113
- Distinguish a Beta-gamma-proteobacteria from an Epsilonproteobacteria
    - B000064, B000076, B000087, B000098, B000106 should be there
    - B000098, B000107, B000111, B000114 should not be there

## *Which markers using?*


- 16S

    ![](images/16s_marker_bias.png "Figure 1 from [Liu et al, 2011](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-S2-S4)")

    Taxonomic profile estimated from 16S rRNA targeted sequencing is biased because of copy number variation

    gold standard for long time

    - sensitive to GC compositional bias
    - really conserved sequences: no way to resolve relationships between closely related organisms

- whole-metagenome shotgun sequences

    ![](images/full_genome_marker_bias.png) "Figure 1 from [Liu et al, 2011](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-S2-S4)")


    classification of whole-metagenome shotgun sequences --> biased estimation because of the variations in genome size or copy number

    Need to focus on single-copy gene families --> more accurate estimates of taxonomic abundance 


- Universally conserved genes
    - Phylogenetic marker genes from representatives of complete bacterial genomes
    - Single copy genes within each genome
    - No subjects to HGT
    - Tools / databases
        - AMPHORA ([Wu & Eisen, 2008](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-10-r151)): 31 proteins for 578 complete bacterial genomes 
        - MLTreeMap ([Stark et al, 2010](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-11-461)): 44 gene families from STRING and RefSeq + 4 function reference families from KEGG and STRING
        - AMPHORA2 ([Wu & Scott, 2012](https://academic.oup.com/bioinformatics/article/28/7/1033/210898)): 31 genes of AMPHORA for bacteria + 104 marker genes for 45 archae
            Recalcitrant to lateral gene transfer
        - MetaPhyler ([Liu et al, 2011](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-S2-S4)): 31 genes of AMPHORA from all complete genomes, NCBI nr protein database and 60 draft genomes + phylogenetic marker genes from Archaea
        - PhyloSift ([Darling et al, 2014](https://peerj.com/articles/243/)): ~800 gene families
            - 37 marker gene families (from 40 families from [Wu et al, 2013](https://doi.org/10.1371%2Fjournal.pone.0077033))
                - Large congruent phylogenetic histories
                - 1% of an average bacterial genomes
            - 16S and 18S ribosomal RNA genes
            - Mitochondrial gene families
            - Eukaryote-specific gene families
            - Viral gene families
- Clade-specific genes
    - strongly conserved within the clade's genomes
    - not possessing substantial local similarity with any sequence outside the clade
    - Tools / databases
        - MetaPhlAn ([Segata et al, 2012](https://www.nature.com/articles/nmeth.2066)): 1,221 species with an average of 231 marker CDSs per species + >115,000 markers at higher taxonomic levels
            - 2,887 genomes in IMG
            - 2 million potential markers
            - 400,141 genes most representative of each taxonomic unit
        - MetaPhlAn2 ([Truong et al, 2015](https://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3589.html)):  ~1 million markers (184 $\pm$ 45 for each bacterial species) from >7,500 species (bacteria, archaea, viruses, fungi, protozoa)

## Methods

![](images/marker_gene_analyis.png)


2 steps:
1. Identify reads that hit to one of the markers --> faster assignment (smaller db thant db with full genomes for all species)
2. Classification

2 methods to taxonomic annotation
- sequence similarity between the read and the marker genes
- phylogenetic information
    - longer to calculation
    - greate accuracy
    

### Sequence similarity approach

Example of MetaPhyler ([Liu et al, 2011](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-S2-S4)):

- Search with BLASTX
- Classification of each sequence individually based on its best reference hit


- $Q$: query sequence
- $G$: gene with best hit
- $b$: BLAST bit score
- $L$: High Scoring Pairs (HSP) length
- Classification
    1. Try to classify $Q$ at genus level: Computing bit score cutoff $b_{cut}$ of gene $G$ using the pre-computed linear regression function
    2. If $b \geq b_{cut}$, assignation of genus of $G$ to $Q$
    3. Else, try to classify at higher taxonomic levels

![](images/metaphyler_classification.jpg) "Figure 5 from [Liu et al, 2011](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-S2-S4)")

MetaPhlAn ([Segata et al, 2012](https://www.nature.com/articles/nmeth.2066), [Truong et al, 2015](https://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3589.html))
- Search with BLAST (MetaPhlAn) or Bowtie2 (MetaPhlAn2) 
- Estimation of each clade's relative abundance: $\frac{\text{total number of reads in each clade}}{\text{nucleotide length of its markers}}$


### Phylogenetic information approach

![](images/amphora.jpg "Example of HMMs with [AMPHORA](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-10-r151)")


Phylotyping of a query sequence

![](images/amphora_classification.jpg)


- Database: for each marker
    - Identification of their protein sequences from representative bacterial genomes
    - Alignment
    - Manually edited and masked: 1 for reliably aligned columns and 0 for ambigous columns
    - Generate local profile HMMs from the seed alignment
- Alignement of the query amino acid sequences onto the trusted and fixed seed alignments
- Trimming of the query alignments
- Insertion into the reference tree using a maximum parsimony method

Phylotyping algorithm
- Identification of immediate ancestor $n_0$ and first internal node $n_1$ with $\geq$ 70% bootstrapping support
- Infer taxonomy of the query with the known descendant leaf nodes of $n_1$ + normalized branch length information


## Limits

- Assumption that the relatively small fraction of the metagenome that is homologous to marker genes represents an accurate sampling of the entire taxonomic distribution of the community
- Not appropriate for taxa that do not contain the markers being explored
- Variation accross markers of the accuracy of annotation based on properties of the marker family
- Efforts to identify phylogenetic clade-specific marker genes and expand the phylogenetic diversity represented in genome sequence databases
- Impossible to know the copy number of a gene for a species with an incomplete genome
- Identification of only a few genes per genome: most reads without classification and composition expressed in terms of relative abundance for all recognize taxa

Alternative: using the overlap of MinHash signatures to estimate the similarity of data sets extremely efficiently, e.g. the overlap between all microbial genomes in GenBank and a metagenomics data set

Tools: Mash, sourmash

For large data sets with hundreds of samples on which performing or interpreting metagenomics assembly is impractical, marker-based approaches are currently the method of choice, especially for environments with a substantial fraction of microbial diversity covered by well-characterized sequenced species.



*Who is There?* Assessing Taxonomic Diversity

.center[![:scale 65%](../img/metagenomics/metagenomics_taxonomy_diversity.jpg)]

.footnote[Figure 2 from [Sharpton, 2014](https://www.frontiersin.org/articles/10.3389/fpls.2014.00209/full)]

???
3 different strategy for taxonomic classification

- sequence similarity based methods, which use the results of a sequence similarity search against a database of a reference set of sequences
- sequence composition based methods, which are based on characteristics of their nucleotide composition (e.g. tetranucleotide usage or codon usage)
- marker-based methods which identify species based on the occurrence of certain specific marker sequences

2 possible outputs