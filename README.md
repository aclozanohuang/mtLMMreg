# mtLMMreg

mtLMMreg enables regularized multitrait multilocus linear mixed model estimation.

Given a trait matrix, a genotype matrix and a relatedness matrix, two regularized estimators are proposed to estimate the fixed effect matrix and the covariance matrices of the genetic and noise components of the linear mixed model:
 - mtlmmlasso implements the mtLMM-L1 estimator which imposes variable selection on the fixed effect parameter matrix
 - mtlmmclust implements the mtLMM-clust estimator which imposes variable selection and traitw-wise clustering on the fixed effect parameter matrix

See regmtlmm/mtlmmlasso.py and regmtlmm/mtlmmclust.py for a detailed description of input, output and options for each method.

Note: mtLMMreg builds on the package liMMBo  (https://github.com/HannahVMeyer/limmbo/) whose code is included in the present repository with a minor modification. See folder `limmbo LICENSE and NOTICE' 


## Install 

Install LIMIX with multi-variate support (v1.0.18):
```bash
pip3 install "limix<2"
```

Download or clone the code provided in the mtLLMreg repository


## Requirements

See require.txt


## Examples

Examples are provided in the notebook demo-mtlmm-reg.
