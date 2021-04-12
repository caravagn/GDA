Genomic Data Analytics (GDA)
================
Giulio Caravagna (<gcaravagna@units.it>)
3/5/2021

*MSc program in Data Science and Scientific Computing. University of
Trieste, Italy*

-   3CFU - 24h, 12 lecture, 2 hours each. 50% theoretical lecture, 50%
    practical session (40’ each).
-   GitHub: <https://github.com/caravagn/GDA>

## Lecturers

-   Dr Giulio Caravagna, [Cancer Data Science
    Laboratory.](http://www.caravagnalab.org)

**Invited guest lecturers**

-   Dr Riccardo Bergamin, University of Trieste
-   Dr Alex Graudenzi, CNR.
-   Dr Salvatore Milite, University of Trieste
-   Dr Daniele Ramazzotti, University of Milan-Bicocca.

------------------------------------------------------------------------

## Program

| Lecturer         | Title                                         | When  |
|:-----------------|:----------------------------------------------|-------|
| Caravagna        | *Variant calling from bulk sequencing*        | 9/4   |
| Caravagna        | *Measuring aneuploidy from bulk sequencing*   | 12/4  |
| Caravagna        | *Integrated quality control of somatic calls* | 16/4  |
| Bergamin         | *Population genetics for cancer*              | 19/4  |
| Caravagna        | *Tumour subclonal deconvolution*              | 21/4  |
| Ramazzotti       | *Somatic mutational signatures*               | 23/4  |
| Bergamin, Milite | *Basics of Single-cell RNA analysis*          | 30/4 |
| Graudenzi        | *Longitudinal evolution from single cell*     | 3/5 |
| Milite           | *Count-based models for single-cell data*     | 5/5 |
| Caravagna        | *Evolutionary based stratifications*          | 10/5|
| Caravagna        | *Population-level models*                     | 12/5 |

------------------------------------------------------------------------

### Preamble

-   Course presentation (https://www.dropbox.com/s/vl3u4gjmitgvmk9/Primer_annotated_lecture.pdf?dl=0):

    -   Cancer Evolution,
    -   Modern Genomics,
    -   Single-cell.
    -   Research at the CDSLab (www.caravagnalab.org).

------------------------------------------------------------------------

### Part 1 - Somatic calling from bulk sequencing

------------------------------------------------------------------------

------------------------------------------------------------------------

**Lecture:** *Variant calling from bulk sequencing* (https://www.dropbox.com/s/2ngbbxiudux8h9v/Somatic_calling_annotated_lecture.pdf?dl=0)

-   (Theory) Mutation calling:

    -   Tumour matched-normal design,
    -   High-level design of GATK
    -   Joint calling model

-   (Practice) Example VCF and PCAWG:

    -   VCF manipulation
    -   27 PCAWG cases (mutation types, burden, etc.)

-   Readings (https://www.dropbox.com/s/bx2kam7tlf7tl5x/readings.zip?dl=0)

    -   (tool) Roth, Andrew, et al. “JointSNVMix: a probabilistic model
        for accurate detection of somatic mutations in normal/tumour
        paired next-generation sequencing data.” Bioinformatics 28.7
        (2012): 907-913.
    -   (tool) Kim, Sangtae, et al. “Strelka2: fast and accurate calling
        of germline and somatic variants.” Nature methods 15.8 (2018):
        591-594.
    -   (tool) Benjamin, David, et al. “Calling somatic snvs and indels
        with mutect2.” BioRxiv (2019): 861054.
    -   (tool) Rimmer, Andy, et al. “Integrating mapping-, assembly-and
        haplotype-based approaches for calling variants in clinical
        sequencing applications.” Nature genetics 46.8 (2014): 912-918.
    -   (tool) GATK (Broad Institute)
       - Training: www.csc.fi/en/web/training/-/gatk2019
       -   Lectures: https://www.youtube.com/watch?v=sM9cQPWwvn4&list=PLjiXAZO27elDHGlQwfd06r7coiFtpPkvI

------------------------------------------------------------------------

**Lecture:** *Measuring aneuploidy from bulk sequencing*

-   (Theory) Aneuploidy and Copy Number calling:

    -   Motivation
    -   ASCAT model
    -   Segmentation

-   (Practice) Example runs with different tools:

    -   ASCAT
    -   Sequenza (inspection of alternative solutions)
    -   Circular Binary Segmentation
    -   Cohort-level distribution of CNAs per chromosome (length,
        percentage, copy state).

-   Readings (https://www.dropbox.com/s/ohz5f7e51dwpg71/readings.zip?dl=0)

    -   (tool) Favero, Francesco, et al. “Sequenza: allele-specific copy
        number and mutation profiles from tumor sequencing data.” Annals
        of Oncology 26.1 (2015): 64-70.
    -   (tool) Van Loo, Peter, et al. “Allele-specific copy number
        analysis of tumors.” PNAS 107.39 (2010): 16910-1691
    -   (tool) Ross, Edith M., et al. “Allele-specific multi-sample copy
        number segmentation in ASCAT.” Bioinformatics (2020).
    -   (tool) Olshen, Adam B., et al. “Circular binary segmentation for
        the analysis of array‐based DNA copy number data.” Biostatistics
        5.4 (2004): 557-572.
    -   (Review) Ben-David, Uri, and Angelika Amon. “Context is
        everything: aneuploidy in cancer.” Nature Reviews Genetics 21.1
        (2020): 44-62
    -   (Review) Weaver, Beth AA, and Don W. Cleveland. “Does aneuploidy
        cause cancer?.” Current opinion in cell biology 18.6 (2006):
        658-667.
    -   (In vivo measurements) Bolhaqueiro, Ana CF, et al. “Ongoing
        chromosomal instability and karyotype evolution in human
        colorectal cancer organoids.” Nature Genetics 51.5 (2019):
        824-834.
    -   (coding) DNAcopy: A Package for Analyzing DNA Copy Data
        <https://bioconductor.org/packages/release/bioc/vignettes/DNAcopy/inst/doc/DNAcopy.pdf>
    -   (coding)Total copy-number segmentation using CBS.
        <https://cran.r-project.org/web/packages/PSCBS/vignettes/CBS.pdf>

------------------------------------------------------------------------

**Lecture:** *Integrated quality control of somatic calls*

-   (Theory) Validating mutations, copy number and tumour purity:

    -   Cancer Cell Fractions
    -   CNAqc
    -   Tumour in Normal contamination (ideas)

-   (Practice) Quality-control of Whole Genome Sequencing data:

    -   Implement mapping of SNVs to CNA segments
    -   Validate MSeq calls with CNAqc
    -   Identify good and bad samples from PCAWG

-   Readings

    -   Househam, Jacob, William CH Cross, and Giulio Caravagna. “A
        fully automated approach for quality control of cancer mutations
        in the era of high-resolution whole genome sequencing.” bioRxiv
        (2021).
    -   Cmero, Marek, et al. “Inferring structural variant cancer cell
        fraction.” Nature communications 11.1 (2020): 1-15.
    -   Yuan, Ke, et al. “Ccube: a fast and robust method for estimating
        cancer cell fractions.” bioRxiv (2018): 484402.

------------------------------------------------------------------------

# Parte 2 - modelling and inference from bulk

------------------------------------------------------------------------

------------------------------------------------------------------------

**Lecture (R Bergamin):** *Population genetics models of growth*

-   (Theory) Branching processes and other models

    - Cancer Evolution as Stochastic Branching Process
    - Markov System and Master equation
    - Some Examples: Moran Model, Wright-Fisher Model, Coalescence
    - Birth-Death Process
    - Luria-Delbruck Model
    - Theory of 1/f tail
    - Quantify Cancer Evolution from VAF Spectrum
    -   Spatial Tumor Growth
    
-   (Practice) Tumour growth simulation:

    -   Simulations of a Branching process and VAF spectrum
        (<https://t-heide.github.io/TEMULATOR/>),
    -   Example tumours from CHESS
        (<https://github.com/sottorivalab/CHESS.cpp>)

-   Readings

    -  TBD

------------------------------------------------------------------------

**Lecture:** *Tumour subclonal deconvolution*

-   (Theory) Subclonal deconvolution:

    -   Tail modelling versus subclones
    -   Read counts analysis
    -   Multi-sample deconvolution

-   (Practice) Deconvolution in practice

    -   MOBSTER and BMix with single-sample data
    -   MOBSTER and VIBER with multi-region data.

-   Readings

    -   Roth, Andrew, et al. “PyClone: statistical inference of clonal
        population structure in cancer.” Nature methods 11.4 (2014):
        396-398.
    -   Gillis, Sierra, and Andrew Roth. “PyClone-VI: scalable inference
        of clonal population structures using whole genome data.” BMC
        bioinformatics 21.1 (2020): 1-16.
    -   Miller, Christopher A., et al. “SciClone: inferring clonal
        architecture and tracking the spatial and temporal patterns of
        tumor evolution.” PLoS Comput Biol 10.8 (2014): e1003665.
    -   Caravagna, Giulio, et al. “Subclonal reconstruction of tumors by
        using machine learning and population genetics.” Nature Genetics
        52.9 (2020): 898-907.
    -   Caravagna, Giulio, et al. “The MOBSTER R package for tumour
        subclonal deconvolution from bulk DNA whole-genome sequencing
        data.” BMC bioinformatics 21.1 (2020): 1-11.

------------------------------------------------------------------------

**Lecture (D Ramazzotti):** *Mutational signatures in human cancers*

- Theory: 
   - Concepts behind mutational signatures 
   - De novo inference of mutational signatures 
   - Solving with non-negative matrix factorization (NMF) 
   - Mutational signature extraction from pan-cancer data 

- Practice: 
   - Examples and best practice on real data 
   - Analysis of breast cancer data 

Readings: 

   - (Concepts on mutational signatures) Alexandrov, Ludmil B., et al. "Signatures of mutational processes in human cancer." Nature 500.7463 (2013): 415-421.
   - (Concepts on mutational signatures) Alexandrov, Ludmil B., et al. "The repertoire of mutational signatures in human cancer." Nature 578.7793 (2020): 94-101.
   - (Tool - SigProfiler) Alexandrov, Ludmil B., et al. "Deciphering signatures of mutational processes operative in human cancer." Cell reports 3.1 (2013): 246-259.
   - (Tool - SparseSignatures) Lal, A., et al. "De Novo Mutational Signature Discovery in Tumor Genomes using SparseSignatures." (2020).
   - (Statistics - Non-negative matrix factorization) Brunet, Jean-Philippe, et al. "Metagenes and molecular pattern discovery using matrix factorization." Proceedings of the national academy of sciences 101.12 (2004): 4164-4169.
   - (Statistics - Non-negative matrix factorization) Owen, Art B., and Patrick O. Perry. "Bi-cross-validation of the SVD and the nonnegative matrix factorization." The annals of applied statistics 3.2 (2009): 564-594.

------------------------------------------------------------------------

# Part 3 - Single-cell genomics transcriptomics

------------------------------------------------------------------------

------------------------------------------------------------------------

**Lecture (R Bergamin, S Milite):** *Basics of Single-cell RNA analysis*

-   (Theory) Single-cell RNA analysis

    -   …
    -   …

-   (Practice) Seurat analysis

    -   …
    -   …

-   Readings

    -   …
    -   …

------------------------------------------------------------------------

Lecture

**Lecture (A Graudenzi):** *Longitudinal evolution from single cell*

-   (Theory) Principles of evolution from single-cell

    -   …
    -   …

-   (Practice) Inference in practice

    -   …
    -   …

-   Readings

    -   …
    -   …


------------------------------------------------------------------------

**Lecture (S Milite):** *Count-based models for single-cell data*

- Theory (approx 75/80 min):
   - Generative modelling as an alternative to pipelines
   - Poisson and Negative binomial distributions
   - Count based modelling, RNA-seq vs scRNA-seq
   - Count models for normalisation (scTransform)
   - Scaling NB models with variational autoencoders (scVI)
   - CONGAS (genotype CNV from scRNA-seq)
   - Elements of Gradient based variational inference
   - Discrete Latent Variable modelling

- Practice (approx 20-30 min):
   - Example run of CONGAS on breast cancer 10x dataset

- Readings:
   - Hafemeister, Christoph, and Rahul Satija. ‘Normalization and Variance Stabilization of Single-Cell RNA-Seq Data Using Regularized Negative Binomial Regression’. Genome Biology 20, no. 1 (23 December 2019): 296. https://doi.org/10.1186/s13059-019-1874-1.
   - Jang, Eric, Shixiang Gu, and Ben Poole. ‘Categorical Reparameterization with Gumbel-Softmax’. ArXiv:1611.01144 [Cs, Stat], 5 August 2017. http://arxiv.org/abs/1611.01144.
   - Kingma, Diederik P., and Max Welling. ‘Auto-Encoding Variational Bayes’. ArXiv:1312.6114 [Cs, Stat], 1 May 2014. http://arxiv.org/abs/1312.6114.
   - Lopez, Romain, Jeffrey Regier, Michael B. Cole, Michael I. Jordan, and Nir Yosef. ‘Deep Generative Modeling for Single-Cell Transcriptomics’. Nature Methods 15, no. 12 (December 2018): 1053–58. https://doi.org/10.1038/s41592-018-0229-2.
   - Milite, Salvatore, Riccardo Bergamin, and Giulio Caravagna. ‘Genotyping Copy Number Alterations from Single-Cell RNA Sequencing’. BioRxiv, 1 January 2021, 2021.02.02.429335. https://doi.org/10.1101/2021.02.02.429335.
   - Sarkar, Abhishek, and Matthew Stephens. ‘Separating Measurement and Expression Models Clarifies Confusion in Single Cell RNA-Seq Analysis’. BioRxiv, 1 January 2020, 2020.04.07.030007. https://doi.org/10.1101/2020.04.07.030007.
   - Schulman, John, Nicolas Heess, Theophane Weber, and Pieter Abbeel. ‘Gradient Estimation Using Stochastic Computation Graphs’. ArXiv:1506.05254 [Cs], 5 January 2016. http://arxiv.org/abs/1506.05254.
   - Svensson, Valentine. ‘Droplet ScRNA-Seq Is Not Zero-Inflated’. Nature Biotechnology 38, no. 2  (February 2020): 147–50. https://doi.org/10.1038/s41587-019-0379-5.
   - 
------------------------------------------------------------------------

------------------------------------------------------------------------

# Part 4 - Population-level inference

------------------------------------------------------------------------

**Lecture:** *Evolutionary based stratifications*

-   (Theory) Detecting repeated evolution from multi-region bulk
    sequencing

    -   Clone-trees and tree expansion
    -   Expectation Maximisation for latent model discovery
    -   Evolutionary distance and cluster

-   (Practice) Inference in practice

    -   Colorectal adenomas with REVOLVER
    -   TRACERx Lung Adencarcinomas with REVOLVER

-   Readings

    -   Caravagna, Giulio, et al. “Detecting repeated cancer evolution
        from multi-region tumor sequencing data.” Nature methods 15.9
        (2018): 707-714.

**Lecture:** *Population-level models*

-   (Theory) Bayesian Networks models

    -   Conjunctive Bayesian Networks
    -   Suppes’ probabilistic causation

-   (Practice) Inference in practice

    -   Analysis of CODREAD with PICNIC
    -   Analysis of other cbio data

-   Readings

    -   Beerenwinkel, Niko, Nicholas Eriksson, and Bernd Sturmfels.
        “Conjunctive bayesian networks.” Bernoulli (2007): 893-909.
    -   Gerstung, Moritz, et al. “Quantifying cancer progression with
        conjunctive Bayesian networks.” Bioinformatics 25.21 (2009):
        2809-2815.
    -   Caravagna, Giulio, et al. “Algorithmic methods to infer the
        evolutionary trajectories in cancer progression.” Proceedings of
        the National Academy of Sciences 113.28 (2016): E4025-E4034.
    -   Ramazzotti, Daniele, et al. “CAPRI: efficient inference of
        cancer progression models from cross-sectional data.”
        Bioinformatics 31.18 (2015): 3016-3026.
    -   Loohuis, Loes Olde, et al. “Inferring tree causal models of
        cancer progression with probability raising.” PloS one 9.10
        (2014): e108358.
