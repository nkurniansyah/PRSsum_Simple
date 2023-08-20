## Introduction

This repository provides information regarding the construction of a
multi weighted and un-weighted polygenic risk score (PRSsum).

Suppose that we have multiple PRSs, all using the same *p* variants\*,
and we sum them with weights. The formula for the *k*<sub>*t**h*</sub>
PRS:

$$\\frac{1}{2p}\\sum\_{i=1}^{p} g\_{i}\\beta\_{ki}$$

Where note that is standardized according to the number of potential
alleles (2 times the number of variants).

For wighted PRSsum, Suppose that we computed the mean and standard
deviation of the PRSs, and they are given by
*μ*<sub>*k*</sub>,*σ*<sub>*k*</sub> for the *K* PRS. Noting that the
number of variants can change between PRSs, we denote by
*p*<sub>*k*</sub> the number of variants in the *k*<sub>*t**h*</sub>
PRS. When summing *K* PRSs with weights, the expression is:

$$wPRSsum=w\_{1}\\frac{\\frac{1}{2p\_{1}}\\sum\_{i=1}^{p}g\_{i}\\beta\_{1i}-\\mu\_{1}}{\\sigma\_{1}}+ \\dots+w\_{1}\\frac{\\frac{1}{2p\_{k}}\\sum\_{i=1}^{p}g\_{k}\\beta\_{ki}-\\mu\_{k}}{\\sigma\_{k}}$$

This can alternatively be written as a single weighted PRS formula after
some rearrangement of terms:
$$wPRSsum=g\_{i}\\left(\\frac{w\_{1}}{2p\_{1}\\sigma\_{1}} \\beta\_{1i} + \\dots + \\frac{w\_{k}}{2p\_{k}\\sigma\_{k}} \\beta\_{ki}\\right)-\\left( \\frac{w\_{1}\\mu\_{1}}{2p\_{1}\\sigma\_{1}}+\\dots+\\frac{w\_{k}\\mu\_{k}}{2p\_{k}\\sigma\_{k}} \\right)$$

For Unweighted PRSsum, we will compute PRSsum without the weight (*w*)

$$PRSsum=\\frac{\\frac{1}{2p\_{1}}\\sum\_{i=1}^{p}g\_{i}\\beta\_{1i}-\\mu\_{1}}{\\sigma\_{1}}+ \\dots+\\frac{\\frac{1}{2p\_{k}}\\sum\_{i=1}^{p}g\_{k}\\beta\_{ki}-\\mu\_{k}}{\\sigma\_{k}}$$

This can alternatively be written as a single unweighted PRS formula
after some rearrangement of terms:
$$PRSsum=g\_{i}\\left(\\frac{1}{2p\_{1}\\sigma\_{1}} \\beta\_{1i} + \\dots + \\frac{1}{2p\_{k}\\sigma\_{k}} \\beta\_{ki}\\right)-\\left( \\frac{\\mu\_{1}}{2p\_{1}\\sigma\_{1}}+\\dots+\\frac{\\mu\_{k}}{2p\_{k}\\sigma\_{k}} \\right)$$

Finnaly, We provide example data set in ./Example directory and code to
generate weighted and unweighted PRSsum.

## Example to construct weighted PRSsum

In ./Example directory, we provide 4 summary statistics from diffrent
GWAS together with scaling and weight files.

In this example: we assiggne AFR as Summary\_stat1, EUR as
Summary\_stat2, FINNGEN as Summary\_stat3 and HIS as Summary\_stat4.

for ./Example/2023-07-17\_PRS\_scaling\_AdjustedBMI\_GWAS.csv file. this
file contain Study name, Mean and SD.

    ##           Study      Mean       SD
    ## 1 Summary_stat1  2.66e-07 1.46e-07
    ## 2 Summary_stat2  6.21e-07 1.95e-07
    ## 3 Summary_stat3 -2.41e-08 4.95e-08
    ## 4 Summary_stat4 -5.62e-07 1.27e-07

for ./Example/2023-07-17\_MGB\_weight\_AdjustedBMI\_GWAS.csv file. this
file contain Study name, prs\_effect.

    ##           Study prs_effect
    ## 1 Summary_stat1 0.04180444
    ## 2 Summary_stat2 0.72934181
    ## 3 Summary_stat3 0.12193262
    ## 4 Summary_stat4 0.04030905

We will use ./Code/create\_wPRSsum.R to construct weighted PRSsum.

    source("Code/create_wPRSsum.R")

    Summary_stat1<-"Example/AFR.txt"
    Summary_stat2<-"Example/EUR.txt"
    Summary_stat3<-"Example/FINNGEN.txt"
    Summary_stat4<-"Example/HIS.txt"

    n_sample_Summary_stat1<- 873501
    n_sample_Summary_stat2<- 938203
    n_sample_Summary_stat3<- 930315
    n_sample_Summary_stat4<- 924899

    scaling_file<-"Example/2023-07-17_PRS_scaling_AdjustedBMI_GWAS.csv"
    weight_file<-"Example/2023-07-17_MGB_weight_AdjustedBMI_GWAS.csv"


    wprssum<- create_wprsum(Summary_stat1=Summary_stat1, Summary_stat2=Summary_stat2,
                            Summary_stat3=Summary_stat3, Summary_stat4=Summary_stat4, 
                            n_sample_Summary_stat1= n_sample_Summary_stat1,
                            n_sample_Summary_stat2=n_sample_Summary_stat2, 
                            n_sample_Summary_stat3=n_sample_Summary_stat3,
                            n_sample_Summary_stat4=n_sample_Summary_stat4,
                            scaling_file=scaling_file, weight_file=weight_file)

    ## Warning in read.table(file = file, header = header, sep = sep, quote
    ## = quote, : incomplete final line found by readTableHeader on 'Example/
    ## 2023-07-17_PRS_scaling_AdjustedBMI_GWAS.csv'

    ## Warning in read.table(file = file, header = header, sep = sep, quote
    ## = quote, : incomplete final line found by readTableHeader on 'Example/
    ## 2023-07-17_MGB_weight_AdjustedBMI_GWAS.csv'

    head(wprssum)

    ##   CHR      rsID    POS A1 A2          BETA
    ## 1   1 rs4040617 843942  A  G -9.577557e-05
    ## 2   1 rs4970383 903175  C  A -6.206930e-05
    ## 3   1 rs4475691 911428  C  T -2.095398e-04
    ## 4   1 rs1806509 918574  C  A  1.875377e-04
    ## 5   1 rs7537756 918870  A  G -2.488695e-04
    ## 6   1 rs1110052 938178  G  T -1.022256e-04
