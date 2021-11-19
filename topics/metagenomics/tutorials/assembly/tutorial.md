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

Reconstruction of 2,631 draft metagenome-assembled genomes from the global oceans

![](images/mag_sampling.jpg "Figure 1 from [Tully et al, 2017](https://www.nature.com/articles/sdata2017203)")

???
TARA ocean data


Reconstruction of 2,631 draft metagenome-assembled genomes from the global oceans

![](images/mag_reconstructed.jpg "Figure 2 from [Tully et al, 2017](https://www.nature.com/articles/sdata2017203)")


![](images/mag_sampling.jpg)

Other applications
- Qualitative understanding of the physiology of the uncultivated microbes
- [Sangwan et al, 2016](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-016-0154-5)
- Examples
    - Identification of the enzymes used for oil and paraffin degradation by Smithella spp ([Tan et al](http://genomea.asm.org/content/2/5/e01085-14.short), [Wawrik et al, 2016](http://onlinelibrary.wiley.com/doi/10.1111/1462-2920.13374/full))
    - insights into metabolic pathways and interactions between microbes in methanogenic bioreactors [Nobu et al 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4511927/)


> **Metagenome-assembled genome (MAG)**: genome assembled using reads from metagenomics datasets

![](images/MetagenomeAssembly.jpg)


# Assembly using De Bruijn graphs

![](images/assembly_de_bruijn.png "Figure from [Vollmers et al, 2017](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169662)")


DBG: one looks for multiple paths through the graph that collectively ‘explain’ all the edges

- Single $k$-mer
    - 2 approaches
        -  MetaVelvet ([Namiki et al, 2012](https://doi.org/10.1093/nar/gks678))
            - Decomposition of the single de Bruijm graph into multiple subgraph (ideally corresponding to different organisms) based on coverage information and graph connectivity
        - Ray Meta ([Boisvert et al, 2012](https://doi.org/10.1186/gb-2012-13-12-r122))
            - Construction contigs by a heuristics-guided graph traversal
    - Importance of k
        - Small k: more sensitive in making connections, but fail to resolve repeats
        - Large k: miss connections and are more sensitive to sequencing errors, but usually create longer contigs

- Iteratively constructed and refined de Bruijn graphs using multiple $k$-mer lengths
    - IDBA
        - Going from small k's to large k's, replacing reads with preassembled contigs at each iteration
    - IDBA-UD ([Peng et al, 2012](https://doi.org/10.1093/bioinformatics/bts174))
        - Tolerate uneven depth of coverage
        - Idea
            1. Generation of a de Bruijn graph from the reads using small k-mers (by default k = 20)
            2. Error correction
            3. Extraction of contigs that are used as ‘reads’ in the graph construction with the next-higher k-mer size
        - Detection of erroneous k-mers and k-mers from different genomes by looking at deviations from the average multiplicity of k-mers in a contig --> more accurately decompose the de Bruijn graph
    - MetaSPAdes ([Nurk et al, 2017](http://genome.cshlp.org/content/27/5/824.short))
        - similar to IDBA with iterative de Bruijn graph refinement, but keeping the complete read information together with preassembled contigs at each step
        - heurestics for graph simplification, filtering and storage
        - use of "strain-contigs" to inform the assembly of high-quality consensus backbone sequences
        - good for assembling libraries sequenced with different technologies (hybrid assembly)
    - Megahit ([Li et al, 2015](https://doi.org/10.1093/bioinformatics/btv033))
        - Uses a succint de Bruijn graph representation (partitions reads using k-mer abundance patterns?) to reduce the memory requirements
        - Only keeps highly reliable k-mers appearing more than once
        - Implementation of a strategy to recover low-depth edges by taking additional k-mers from high-quality reads, which increases the contiguity of low-depth regions

Several factors impact the performance of DBG assemblers: 
1. sequencing errors: create ‘false’ k-mers
    Every error impacts at most k different k-mers; thus, the impact of sequencing errors increases with the size of k.
    Initial de Bruijn assemblers used spectral correction: use a fixed k-mer count threshold to define ‘correct’ k-mers, strategy that is insufficient in metagenomic data sets with varying coverage levels. Recent approaches to correction have been proposed, which can correct data without assuming uniform coverage
2. repeats: create additional edges in the graph, increasing the number of possible traversals
    The longer the size of k, the fewer nodes in the graph are repetitive
3. the presence of strain variants: imilar challenge as sequencing errors
4. the depth of sequencing coverage: impacts the connectivity of the assembly graph
    A path stretching from a read to the next until it covers an entire genome can only be found if adjacent reads share k-mers. At low depths of coverage, the adjacent reads are only expected to overlap by a small extent, and as a result, the assembly is only possible for small values of k.

large values of k
- reduce the complexity of the graph and impact of repeats
- requires longer sequences (longer than the k-mer size) and higher depth of coverage
- leads to an increased impact of sequencing errors

# Overlap Layout Consensus (OLC)

- overlap-layout-consensus graph approach
    - Omega

![](images/assembly_olc.png "Figure from [Vollmers et al, 2017](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169662)")

![](images/DISCO.png "[DISCO](https://disco.omicsbio.org/)")


# The reality

![](images/metagenome_assembly_reality.jpg "Figure from [Iverson et al, 2012](hhttp://science.sciencemag.org/content/335/6068/587.full)")


# Assembly + Binning

![](images/assembly_binning.png "https://www.slideshare.net/dparks1134/parks-kmer-metagenomics") 