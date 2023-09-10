## Introduction

This repository provides information regarding the construction of a
multi weighted and un-weighted polygenic risk score (PRSsum).

Suppose that we have multiple PRSs, all using the same *p* variants\*,
and we sum them with weights. The formula for the *k*<sub>*t**h*</sub>
PRS:

$$\frac{1}{2p}\sum\_{i=1}^{p} g\_{i}\beta\_{ki}$$

Where note that is standardized according to the number of potential
alleles (2 times the number of variants).

For weighted PRSsum, Suppose that we computed the mean and standard
deviation of the PRSs, and they are given by
*μ*<sub>*k*</sub>,*σ*<sub>*k*</sub> for the *k* PRS. Noting that the
number of variants can change between PRSs, we denote by
*p*<sub>*k*</sub> the number of variants in the *k*<sub>*t**h*</sub>
PRS. When summing *K* PRSs with weights, the expression is:

$$wPRSsum=w\_{1}\frac{\frac{1}{2p\_{1}}\sum\_{i=1}^{p\_1} g\_{i}\beta\_{1i}-\mu\_{1}}{\sigma\_{1}}+ \dots+w\_{1}\frac{\frac{1}{2p\_{k}}\sum\_{i=1}^{p\_k}g\_{k}\beta\_{ki}-\mu\_{k}}{\sigma\_{k}}$$

Denote by *p* the total number of variants used in at least one of the
component PRSs, where the corresponding variant weight *β* may equal
zero in a component PRS that does not include this variant. The above
expression can alternatively be written as a single weighted PRS formula
after some rearrangement of terms:
$$wPRSsum=\sum\_{i=1}^p g\_{i}\left(\frac{w\_{1}}{2p\_{1}\sigma\_{1}} \beta\_{1i} + \dots + \frac{w\_{k}}{2p\_{k}\sigma\_{k}} \beta\_{ki}\right)-\left( \frac{w\_{1}\mu\_{1}}{2p\_{1}\sigma\_{1}}+\dots+\frac{w\_{k}\mu\_{k}}{2p\_{k}\sigma\_{k}} \right)$$

For Unweighted PRSsum, we will compute PRSsum without the weight (*w*)

$$PRSsum=\frac{\frac{1}{2p\_{1}}\sum\_{i=1}^{p\_1}g\_{i}\beta\_{1i}-\mu\_{1}}{\sigma\_{1}}+ \dots+\frac{\frac{1}{2p\_{k}}\sum\_{i=1}^{p\_k}g\_{k}\beta\_{ki}-\mu\_{k}}{\sigma\_{k}}$$

As before, this can alternatively be written as a single unweighted PRS
formula after some rearrangement of terms:
$$PRSsum=\sum\_{i=1}^pg\_{i}\left(\frac{1}{2p\_{1}\sigma\_{1}} \beta\_{1i} + \dots + \frac{1}{2p\_{k}\sigma\_{k}} \beta\_{ki}\right)-\left( \frac{\mu\_{1}}{2p\_{1}\sigma\_{1}}+\dots+\frac{\mu\_{k}}{2p\_{k}\sigma\_{k}} \right)$$

We provide an example data set in ./Example directory and code to
generate weighted and unweighted PRSsum.

## Example to construct weighted PRSsum

Here is an example of how to construct weighted summations of variant
weights from a few PRS weight files.

## 1. Prepare the variants weight files.

First, we must align alleles across the study and consolidate them into
a unified data frame. Below is an illustrative example of how to
accomplish this task. We provided a function for allele matching across
the study, which can be found in the file “./Code/match\_allele.R”.
NOTE: Please ensure that the column names in your dataset match those
used in our provided example to prevent any errors.

    library(data.table)
    library(dplyr)
    library(tidyverse)

    source("Code/match_allele.R")

    #Open all the files and modify the column names to align with the function's requirements.
    #1. Create the list for all variants weight
    studies<-c("AFR","EUR","HIS","FINNGEN")
    study_list<-list()
    for(study in studies){
      variant_weight_file<- paste0("Example/",study,".txt")
      
      variant_weight<-fread(variant_weight_file, data.table = F)
      colnames(variant_weight)<-c("Chromosome", "rsID","Position", "Allele1","Allele2", "BETA")

      study_list[[study]]<- variant_weight
    }

    #example of variants weight
    head(study_list[["AFR"]])

    ##   Chromosome      rsID Position Allele1 Allele2          BETA
    ## 1          1 rs4040617   843942       A       G -1.417655e-05
    ## 2          1 rs4970383   903175       C       A -7.481873e-05
    ## 3          1 rs4475691   911428       C       T -2.898467e-05
    ## 4          1 rs1806509   918574       C       A  2.350129e-04
    ## 5          1 rs7537756   918870       A       G -2.977790e-04
    ## 6          1 rs1110052   938178       G       T  1.247273e-05

Ensure the PRSsum scaling file is in the correct format

    ### PRSsum Scaling
    PRSsum_Scaling<- fread("Example/PRSsum_Scaling.csv", data.table = F)
    PRSsum_Scaling

    ##     Study      Mean       SD N_variants
    ## 1     AFR  2.66e-07 1.46e-07     800000
    ## 2     EUR  6.21e-07 1.95e-07     810000
    ## 3 FINNGEN -2.41e-08 4.95e-08     820000
    ## 4     HIS -5.62e-07 1.27e-07     850000

Next, we can ensure that all the alleles match across the study by
providing reference SNPs. In this example, we used EUR as the reference
SNP. Below example of the reference SNP files:

    ref_snp<- fread("Example/Reference_snp.txt", data.table = F)
    head(ref_snp)

    ##   Chromosome       rsID Position Allele1_ref Allele2_ref
    ## 1          1  rs4040617   843942           A           G
    ## 2          1  rs4970383   903175           C           A
    ## 3          1  rs4475691   911428           C           T
    ## 4          1  rs1806509   918574           C           A
    ## 5          1  rs7537756   918870           A           G
    ## 6          1 rs28576697   935265           T           C

Next we can match all the variants weight with the reference SNP

    variants_weight_clean<- match_allele(refrence_snp=ref_snp, 
                                         list_variants_weight=study_list,
                                         match_by_position=FALSE)

    ## variant weight will be matched using rsID, Chromosome and Position

    ##  All the allele from AFR are matched with the reference SNP

    ## There are 12 SNP's in AFR are not found in the reference SNPs, These SNPs will be assigned as new reference SNPs

    ## variant weight will be matched using rsID, Chromosome and Position

    ##  All the allele from EUR are matched with the reference SNP

    ## variant weight will be matched using rsID, Chromosome and Position

    ##  All the allele from HIS are matched with the reference SNP

    ## There are 6 SNP's in HIS are not found in the reference SNPs, These SNPs will be assigned as new reference SNPs

    ## variant weight will be matched using rsID, Chromosome and Position

    ##  All the allele from FINNGEN are matched with the reference SNP

    ## There are 2 SNP's in FINNGEN are not found in the reference SNPs, These SNPs will be assigned as new reference SNPs

    head(variants_weight_clean)

    ##   Chromosome Position Allele1 Allele2      rsID           AFR           EUR
    ## 1          1   843942       A       G rs4040617 -1.417655e-05 -2.457128e-05
    ## 2          1   903175       C       A rs4970383 -7.481873e-05  3.581781e-05
    ## 3          1   911428       C       T rs4475691 -2.898467e-05 -1.396116e-05
    ## 4          1   918574       C       A rs1806509  2.350129e-04 -1.237048e-05
    ## 5          1   918870       A       G rs7537756 -2.977790e-04 -2.649647e-05
    ## 6          1   938178       G       T rs1110052  1.247273e-05 -7.567884e-05
    ##             HIS       FINNGEN
    ## 1 -1.689371e-04 -1.082757e-05
    ## 2  6.725460e-06 -9.154979e-05
    ## 3 -5.658432e-05 -1.254618e-04
    ## 4 -3.133284e-06  1.324633e-04
    ## 5 -1.369418e-04 -9.260439e-05
    ## 6  1.190413e-05  3.451193e-05

## 2. Create summation variant weight

If you intend to run weighted variant weights (weighted PRSsum), you
need to provide the file paths for the weight file. Follow the example
below to format the column name in the weight file.

    weight<- fread("Example/Weight.csv", data.table = F)
    weight

    ##     Study prs_effect
    ## 1     AFR 0.04180444
    ## 2     EUR 0.72934181
    ## 3 FINNGEN 0.12193262
    ## 4     HIS 0.04030905

Weighted variant weights (weighted PRSsum) can be generated using the
following command:

    source("Code/create_PRSsum.R")
    weigted_variants_weight<- create_prsum(variant_weights=variants_weight_clean, 
                                           PRSsum_scaling=PRSsum_Scaling, 
                                           weight_file="Example/Weight.csv",
                                           chr_col_name="Chromosome", 
                                           pos_col_name="Position", 
                                           effect_allele_col_name="Allele1",
                                           other_allele_col_name="Allele2",
                                           rsID_col_name="rsID")

    ## Run weighted PRSsum

    head(weigted_variants_weight)

    ##        rsID chr_name chr_position effect_allele other_allele effect_weight
    ## 1 rs4040617        1       843942             A            G -1.070705e-04
    ## 2 rs4970383        1       903175             C            A -6.694654e-05
    ## 3 rs4475691        1       911428             C            T -2.364286e-04
    ## 4 rs1806509        1       918574             C            A  2.118719e-04
    ## 5 rs7537756        1       918870             A            G -2.791236e-04
    ## 6 rs1110052        1       938178             G            T -1.184337e-04

If you intend to construct an un-weighted PRS summation (unweighted
PRSsum), it can be generated using the following command:

    source("Code/create_PRSsum.R")
    unweigted_variants_weight<- create_prsum(variant_weights=variants_weight_clean, 
                                           PRSsum_scaling=PRSsum_Scaling, 
                                           weight_file=NULL,
                                           chr_col_name="Chromosome", 
                                           pos_col_name="Position", 
                                           effect_allele_col_name="Allele1",
                                           other_allele_col_name="Allele2",
                                           rsID_col_name="rsID")

    ## Run unweighted PRSsum

    head(unweigted_variants_weight)

    ##        rsID chr_name chr_position effect_allele other_allele effect_weight
    ## 1 rs4040617        1       843942             A            G -0.0010543249
    ## 2 rs4970383        1       903175             C            A -0.0013034897
    ## 3 rs4475691        1       911428             C            T -0.0019758356
    ## 4 rs1806509        1       918574             C            A  0.0025840996
    ## 5 rs7537756        1       918870             A            G -0.0031336276
    ## 6 rs1110052        1       938178             G            T  0.0002940934
