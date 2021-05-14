GDA Exam 2021
================
Giulio Caravagna (<gcaravagna@units.it>)
5/14/2021

# Group final project

You will carry out a multi-step analysis of **bulk WGS** data, using the
tools that we have seen during the lectures. The aim is to determine
subgroups of patients (i.e., clusters) of similar evolutionary
trajectories, if they exist in the data. The ones that are most
important are prognostic (i.e., associate with survival analysis).

You need to implement the steps discussed below, assuming to work with
data from one PCAWG tumour type (e.g., colon). Ideally, you should
analyse multiple, if not all, tumour types.

### What data you will work with

You will be given an extended version of PCAWG data, with one `rds` file
per patient formatted in the usual way. The cohort metadata is available
as a [CSV
file](https://www.dropbox.com/s/97fqh6lt0g8qf5f/pcawg_icgc_cohort.csv?dl=0).

In these data we provide you with *somatic driver annotations* derived
from simple somatic mutations (SNVs); this type of information was not
available in the batch we shared during lectures. Drivers are annotated
as two extra columns in the data matrix, `is_driver` and `driver_label`,
which are `TRUE` if the event is a putative driver.

In one step of the analysis, you will have to augment a little bit this
drivers list.

### Quality control somatic calls

The aim is to identify *reasonable criteria* to determine what samples
should be forwarded for downstream analysis, using the available QC
scores or any other (automated) criterion you want to create from a
CNAqc analysis.

For instance, you could use hard cutoffs on CNAqc analysis (e.g., take
the QC statuses by the tool), or perform some more complicated decision
task using some form of classifier. A possibility is a logistic
classifier, which you would supervise annotating `PASS`/`FAIL` status by
manual assessment of some hundred cases, that determines cuts for
inclusion based on covariates such as tumour purity, average coverage,
coverage per karyptype, mutations per karyotype, etc.

For this part of the analysis, you should prepare to show some `PASS`
and `FAIL` cases to motivate your criteria.

### Perform tumour subclonal deconvolution

We want to perform a CCF-based deconvolution using MOBSTER, and
determine which clusters will be used for the next analysis.

CCF values can be computed from CNAqc, following the steps of [this
vignette](https://caravagnalab.github.io/CNAqc/articles/a4_ccf_computation.html).
You should take a decision about what karyotypes you will use, based on
the previous QC step (e.g., a `FAIL` karyotype might be removed). The
CNAqc CCF computation will also return a `PASS` or `FAIL` status per
karyotype, you should take a decision also based on this QC-metric.

Once you have selected CCF values that you consider reliable, you can
fit them with MOBSTER. To do that, you need to add a `VAF` column to
your data - there is one already - and put the CCF values in there. The
deconvolution can select only mutations with VAF above a cutoff (e.g.,
`0.05`), and set `auto_setup = "FAST"` in
[mobster\_fit](https://caravagnalab.github.io/mobster/reference/mobster_fit.html)
to get faster runs.

Investigate critically the quality of the fits obtained by your
analysis. From every fit, you should flag tail mutations and determine
the growth parameters of the tumour. Tail mutations can be determined
from clustering assignments using
[Clusters](https://caravagnalab.github.io/mobster/reference/Clusters.html)
function from MOBSTER. The evolutionary parameters are computed by
function
[evolutionary\_parameters](https://caravagnalab.github.io/mobster/reference/evolutionary_parameters.html),
also in MOBSTER. Measure the mutation rate *μ* for the whole biopsy and
the selection coefficient *s* or the age *t* for each Beta subclone in
MOBSTER - recall that `C2`, `C3`, …, are subclonal, while `C1` is
clonal.

### Cluster evolutionary trajectories

Here we want to compute evolutionary subgroups in the cancer type
selected, and determine if they have significant survival probability
using PCAWG clinical data (where available).

From non-tail mutations, you can prepare [REVOLVER
inputs](https://caravagnalab.github.io/revolver/articles/Input_formats.html)\_.
In this case call the biopsy `Primary` to format the `CCF` string, e.g.,
`Primary:.54` is a 54% VAF in the sample. For each mutation, use the CCF
of its MOBSTER cluster (Beta mean), instead of the mutation’s own CCF.

To run REVOLVER we need to consolidate driver annotations. In
particular, we want to distinguish between *tumour suppressor genes*,
and *oncogenes*. A list of the roles for each possible gene is
accessible as Tier-1 data from the [Cancer Gene Census
(CGC)](https://cancer.sanger.ac.uk/census), scrolling to the end of the
page and downloading the TSV table. Note that some of the drivers you
might find already annotated in the PCAWG data are not present in the
CGC, or not have a suppressor/oncogene status associated. For those, you
should do nothing (just keep them as they are).

For suppressors, instead, you want to check if the gene, besides the
mutation present in the data, harbours *also* an LOH in the segment
where the gene maps. If it does, the `driver_label` in the data should
distinguish suppressors with just a mutation (e.g., `TP53 mut`), from
those with complete inactivation (mutation plus LOH, e.g.,
`TP53 inactive`). In order to map genes to copy number segment you can
use geen coordinates available in [this
tibble](https://github.com/Militeee/rcongas/blob/master/data/hg19_gene_coordinates.rda).

For oncogenes, instead, you want to check if the gene has mutations and
*sits on top* of an amplified copy number segment (`2:1` or `2:2`). This
information is implicitly available from CCF computation: every oncogene
with mutation multiplicity *m* &gt; 1 is on the amplified chromosome
segment. Split your oncogenes as just mutated (e.g., `PIK3CA mut`), or
mutation and amplified, e.g., `PIK3CA mutamp`.

An important result of this selection step is that you will be losing
many drivers if they map in segments that are not supported by CNAqc.
Try to keep track of that information and reason critically on the part
of signal that you are missing by considering only karyotypes `1:0`,
`1:1`, `2:0`, `2:1` and `2:2`.

You can run REVOLVER with default parameters, following the [TRACERx
vignette](https://github.com/caravagn/revolver.misc/blob/master/vignette_TRACERx_Hanjani_et_al/vignette_TRACERx_Hanjani_et_al.md).
You will apply cohort-level frequency filters set up to retain enough
driver events (% in the cohort), which you need to correlate
trajectories with REVOLVER. Exclude the bootstrap from your analysis, or
leave that as an option to compute at a later stage.

### Perform surival analysis

Here we want to test if some of the groups show distinct trends of
survival, using standard methods (Kaplan-Meier curves etc.) available in
R.

You will perform survival analysis between pairs of REVOLVER clusters,
and across all clusters. To do that, you will need [clinical
data](https://dcc.icgc.org/api/v1/download?fn=/PCAWG/clinical_and_histology/pcawg_donor_clinical_August2016_v9.xlsx)
from PCAWG.

You might miss some patients where the data are not available in the
2018 release batch; try to see if the more updated
[repository](https://dcc.icgc.org/) has more recent data. To query the
web portal, you will need to convert the 16 digits sample id to the ICGC
one, using columns `icgc_donor_id` and `aliquot_id` in the cohort
metadata table. In general, you might have that for some groups you have
little surival data to use.

Survival analysis can be done following [this
vignette](https://github.com/caravagn/revolver.misc/blob/master/vignette_survival_Breast_Yates_et_al/vignette_survival_Breast_Yates_et_al.md).
Use `donor_survival_time` (number) and `donor_vital_status`
(categorical) to predict survival. A key difference from that vignette
and this study, is that here you have surival data for the same patients
you have REVOLVER clusters. In the other no, and an external dataset was
used to classify patients based on the clusters’ enriched drivers. You
can therefore avoid that part in your case study.

# Assessment

The project can be carried out in groups of one or more students. Each
of you must have a distinct role defined from the very beginning.

You can divide the tasks among the team, or participate all together in
every task (srrange this the way you prefer). You will be mentored for
all the project (weekly), where each one of you will have to present his
progress and discuss doubts and technical issues, where they raised.

You should work like a team, organise a GitHub repository and create
reports of your results (images, tables with statistics etc.). In the
end, you will be giving a presentation, with slides, about your most
salient results. Try to look critically to your data and the analysis
you are carrying out.
