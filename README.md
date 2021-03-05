Genomic Data Analytics (GDA)
================
Giulio Caravagna (<gcaravagna@units.it>)
3/5/2021

*MSc program in Data Science and Scientific Computing. University of
Trieste, Italy*

-   3CFU - 24h, 12 lecture, 2 hours each. 50% theoretical lecture, 50%
    practical session (40’ each).
-   GitHub: <https://github.com/caravagn/GDA>

### Invited lecturers

-   Dr Alex Graudenzi, CNR.
-   Dr Daniele Ramazzotti, University of Milan-Bicocca.
-   Dr Salvatore Milite, University of Trieste
-   Dr Riccardo Bergamin, University of Trieste

------------------------------------------------------------------------

### Part 1 - Somatic calling from bulk sequencing

------------------------------------------------------------------------

------------------------------------------------------------------------

**Lecture:** *Variant calling from bulk sequencing*

-   (Theory) Introduction to the course:

    -   Cancer Evolution,
    -   Modern Genomics,
    -   Single-cell.
    -   Research at the CDSLab (www.caravagnalab.org).

-   (Theory) Mutation calling:

    -   tumour-matched-normal design,
    -   panel of normals (Binomial testing)
    -   Joint Binomial probabilistic model

-   (Practice) R refresher, example VCF and PCAWG:

    -   Tidyverse
    -   VCF manipulation
    -   27 PCAWG cases (mutation types, burden, etc.)

-   Readings

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
        -   Training: www.csc.fi/en/web/training/-/gatk2019
        -   Lectures:
            <https://www.youtube.com/watch?v=sM9cQPWwvn4&list=PLjiXAZO27elDHGlQwfd06r7coiFtpPkvI>

------------------------------------------------------------------------

**Lecture:** *Measuring aneuploidy from bulk sequencing*

-   (Theory) Aneuploidy and Copy Number calling:

    -   Why it matters
    -   ASCAT model
    -   Segmentation

-   (Practice) Example runs with different tools:

    -   ASCAT
    -   Sequenza (inspection of alternative solutions)
    -   Circular Binary Segmentation
    -   Cohort-level distribution of CNAs per chromosome (length,
        percentage, copy state).

-   Readings

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

# Parte 2 - mathematical modelling e inference from bulk

------------------------------------------------------------------------

------------------------------------------------------------------------

**Lecture (R Bergamin):** *Population genetics models of growth*

-   (Theory) Branching processes and other models

    -   ….
    -   ….
    -   ….

-   (Practice) Tumour growth simulation:

    -   Synthetic tumour generation with TEMULATOR
        (<https://t-heide.github.io/TEMULATOR/>),
    -   Inspecting VAF distributions for subclones that are about to
        sweep, or too small to detect
    -   Example tumours from CHESS
        (<https://github.com/sottorivalab/CHESS.cpp>)

-   Readings

    -   …
    -   …

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

Lecture

    - (Theory) Phylogenetic inference {Neighbour joining, clones trees, mutation trees, SCITE}
    - Sample trees for MSeq
    - Full subclonal deconvolution with multi-region data

------------------------------------------------------------------------

**Lecture (D Ramazzotti):** *Mutational signatures and processes*

-   (Theory) Deconvolution of signatures from SNVs :

    -   …
    -   …

-   (Practice) Deconvolution in practice

    -   …
    -   …

-   Readings

    -   …
    -   …

------------------------------------------------------------------------

# Part 3 - Single-cell genomics transcriptomics

------------------------------------------------------------------------

------------------------------------------------------------------------

Lecture

    - (Theory) Single-cell RNA analysis (Riccardo/Salvatore) {Seurat,Scanpy}
    - Tutorial vignettes from the tools for scRNAseq

Lecture

    - (Theory) Longitudinal evolution single cell (Alex) {SCITE longitudinal, LACE}
    - Longitudinal evolution from scRNAseq

Lecture

    - (Theory) Count-based modelling single-cell RNA (Salvatore) {clone align, CONGAS}
    - CONGAS vignette from scRNAseq

—————

# Part 4 - Population-level inference

Lecture

    - (Theory) Repeated evolution {REVOLVER}
    - Analysis of TRACERx data and MSeq adenocarcinomas

Lecture

    - (Theory) Population-level models {Bayesian Networks models}
    - Analysis of CODREAD with PICNIC
    - Analysis of other cbio/TCGA datasets

—————
