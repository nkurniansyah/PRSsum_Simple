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
