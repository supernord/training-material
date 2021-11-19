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


Questions

- How many species are in a sample?
- How similar are these samples?]

![](images/diversity.png)

???

Ecological and biogeographic measures ([Oulas et al 2015](http://journals.sagepub.com/doi/10.4137/BBI.S12462), [Chiarucci et al 2011](http://rstb.royalsocietypublishing.org/content/366/1576/2426.short))
- species diversity
- species richness
- uniformity of the communities


## $\alpha$ diversity

> **$\alpha$ diversity**: measure of the diversity in a single sample

![](images/alpha_diversity.png)

Richness: quantifies how many different types the dataset of interest contains

![](images/read_abundance_species_abundance.png)

Low correlation between read abundance and species abundance

.footnote[[Edgar, 2017](https://www.biorxiv.org/content/early/2017/04/04/124149)]


> **Rarefaction**: indicator if enough observations have been made to get a good measurement of an $\alpha$ diversity metric

![](images/rarefaction_curve.png)

.footnote[[USEARCH documentation](http://drive5.com/usearch/manual/rare.html)]

???
$R$: measure

rarefaction curve: plot of the value of a measured quantity against the number of observations used in the calculation

Values of R for smaller numbers of observations: random subsets

No convergence: we cannot make a good estimate of R for the full population

2 possibilities
- need more samples to get a good estimate (not yet observed all the taxa present)
- spurious OTUs due to sequencing error (indefinite increase with the number of reads)

---
### Taxonomy and diversity analysis: $\beta$ diversity

> **$\beta$ diversity**: measure of the pairwise dissimilarity of samples

UniFrac metric: phylogenetic-based beta diversity

![](images/unifrac.png)