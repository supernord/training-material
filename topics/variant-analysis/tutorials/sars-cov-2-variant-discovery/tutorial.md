---
layout: tutorial_hands_on

title: Identification of allelic variants in SARS-CoV-2 from deep sequencing reads
level: Introductory
enable: false
zenodo_link: ''
questions:
- Which biological questions are addressed by the tutorial?
- Which bioinformatics techniques are important to know for this type of data?
objectives:
- The learning objectives are the goals of the tutorial
- They will be informed by your audience and will communicate to them and to yourself
  what you should focus on during the course
- They are single sentences describing what a learner should be able to do once they
  have completed the tutorial
- You can use Bloom's Taxonomy to write effective learning objectives
time_estimation: 3H
key_points:
- The take-home messages
- They will appear at the end of the tutorial
contributors:
- wm75
- bebatut

---


# Introduction
{:.no_toc}

Effectively monitoring global infectious disease crises, such as the COVID-19 pandemic, requirescapacity to generate and analyze large volumes of sequencing data in near real time. These data have proven essential for monitoring the emergence and spread of new variants, and for understanding theevolutionary dynamics of the virus.

Two sequencing platforms in combination with several established library preparation strategies are epredominantly used to generate SARS-CoV-2 sequence data. However, data alone do not equal knowledge: they need to be analyzed. The Galaxy community developed analysis workflows to support the **identification of allelic variants (AVs) in SARS-CoV-2 from deep sequencing reads**. 

In this tutorial we will see how to run these workflows for the different type of input data:

- Single end data derived from Illumina-based RNAseq experiments
- Paired end data derived from Illumina-based RNAseq experiments
- Paired-end data generated with Illumina-based Ampliconic (ARTIC) protocols
- ONT fastq files generated with Oxford nanopore (ONT)-based Ampliconic (ARTIC) protocols

> ### Agenda
>
> In this tutorial, we will cover:
>
> 1. TOC
> {:toc}
>
{: .agenda}

# Prepare Galaxy, data and workflows

Any analysis should get their own Galaxy history. So let's start by creating a new one:

> ### {% icon hands_on %} Hands-on: Prepare the Galaxy history
>
> 1. Create a new history for this analysis
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Rename the history
>
>    {% snippet faqs/galaxy/histories_rename.md %}
>
{: .hands_on}

## Import auxiliary datasets

For extracting the variants, we need to get the SARS-CoV-2 reference and gene name aliases (to give gene regions easily recognizable names). In addition, we need for ARTIC protocols 2 extra datasets.


> ### {% icon hands_on %} Hands-on: Import auxiliary datasets
>
> 2. Import the auxiliary datasets from the shared data library / shared history
>
>    - For RNAseq experiments
>    - For ARTIC
>
>    ```
>    
>    ```
>    ***TODO***: *Add the files by the ones on Zenodo here (if not added)*
>
>    ***TODO***: *Remove the useless files (if added)*
>
>    {% snippet faqs/galaxy/datasets_import_via_link.md %}
>
>    {% snippet faqs/galaxy/datasets_import_from_data_library.md %}
>
>
>    <!-- Add tip to say that on usegalaxy.* up-to-date histories are available -->
>
{: .hands_on}

## Get data




## Get the workflows


# From FASTQ to annotated allelic variants



# Report generation



# Conclusion
{:.no_toc}

Sum up the tutorial and the key takeaways here. We encourage adding an overview image of the
pipeline used.