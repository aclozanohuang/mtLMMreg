# mtLMMreg

mtLMMreg enables regularized multitrait multilocus linear mixed model estimation.

Given a trait matrix, a genotype matrix and a relatedness matrix, two regularized estimators are proposed to estimate the fixed effect matrix and the covariance matrices of the genetic and noise components of the linear mixed model:
 - `regmtlmm/mtlmmlasso.py` implements the mtLMM-L1 estimator which imposes variable selection on the fixed effect parameter matrix. See Example 1 in [1]
 - `regmtlmm/mtlmmclust.py` implements the mtLMM-clust estimator which imposes variable selection and trait-wise clustering on the fixed effect parameter matrix. See Example 2 in [1]

Note: mtLMMreg makes use of the package LiMMBo  (https://github.com/HannahVMeyer/limmbo/) whose code is included in the present repository with a minor modification. See folder `limmbo LICENSE and NOTICE` for details. 

[1] "Regularized Multi-trait Multi-locus Linear Mixed Models for Genome-wide Association Studies and Genomic Selection in Crops", Aurelie Lozano, Hantian Ding, Naoki Abe and Alexander E. Lipka, unpublished manuscript, 2022.

## Install 

- Install LIMIX with multi-variate support (v1.0.18):
```bash
pip3 install "limix<2"
```

- Download or clone the code provided in the mtLLMreg repository. Note: Conda and pip installation will be made available once the repository becomes public.


## Requirements

See requires.txt


## Usage
-  `beta,Cg,Cn = mtlmmlasso(filename_pheno, filename_geno, filename_relatedness, maxiter, alpha)`
   
    Arguments:
    
        - filename_pheno (string):       
            name of the csv file containing the phenotype matrix 
            see limmbo/io/test/data/pheno.csv for required format
        - filename_geno (string):        
            name of the csv file containing the genotype matrix 
            see limmbo/io/test/data/genotypes.csv for required format
        - filename_relatedness (string): 
            name of the csv file containing the relatedness matrix 
            see limmbo/io/test/data/relatedness.csv for required format 
        - maxiter (int):
            maximum number of algorithm iterations 
        - alpha (double):    
            L1 penalty parameter controlling the sparsity of the fixed effects matrix (the larger the sparser) 
        
    Returns:
        (tuple):
            tuple containing:
            
            - **beta** (numpy.array):
              fixed effects matrix
            - **Cg** (numpy.array):
              genetic variance component
            - **Cn** (numpy.array):
              noise variance component


-  `beta,Cg,Cn = mtlmmclust(filename_pheno, filename_geno, filename_relatedness, maxiter, C1, C2)`
   
    Arguments:
    
        - filename_pheno (string):       
            name of the csv file containing the phenotype matrix 
            see limmbo/io/test/data/pheno.csv for required format
        - filename_geno (string):        
            name of the csv file containing the genotype matrix 
            see limmbo/io/test/data/genotypes.csv for required format
        - filename_relatedness (string): 
            name of the csv file containing the relatedness matrix 
            see limmbo/io/test/data/relatedness.csv for required format 
        - maxiter (int):
            maximum number of algorithm iterations 
        - C1 (double):    
            L1 penalty parameter controlling the sparsity of the fixed effects matrix  
        - C2 (double):    
            clustering penalty parameter controlling the trait-wise clustering 
            
        
    Returns:
        (tuple):
            tuple containing:
            
            - **beta** (numpy.array):
              fixed effects matrix
            - **Cg** (numpy.array):
              genetic variance component
            - **Cn** (numpy.array):
              noise variance component
          
   
- Examples are provided in the notebook `demo-mtlmm-reg.ipynb`.

