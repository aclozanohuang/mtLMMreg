import pandas as pd
from scipy import linalg
import torch
from sklearn import linear_model
import numpy as np
from sklearn.metrics import mean_squared_error

from pkg_resources import resource_filename

from limmbo.io.reader import ReadData
from limmbo.io.input import InputData
from limmbo.core.vdsimple import vd_reml  

from regmtlmm.twclust import twclust
from numpy import linalg as LA


def mtlmmclust(filename_pheno=None, filename_geno=None, filename_relatedness=None, maxiter=10, C1=1000, C2=20):
    r"""
    Compute the mtlmm-clust estimator

    Arguments:
        filename_pheno (string):       
            name of the csv file containing the phenotype matrix 
            see limmbo/io/test/data/pheno.csv for required format
        filename_geno (string):        
            name of the csv file containing the genotype matrix 
            see limmbo/io/test/data/genotypes.csv for required format 
        filename_relatedness (string): 
            name of the csv file containing the relatedness matrix 
            see limmbo/io/test/data/relatedness.csv for required format
        maxiter (int):
            maximum number of algorithm iterations 
        C1 (double):    
            L1 penalty parameter controlling the sparsity of the fixed effects matrix  
        C2 (double):    
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
            
    """
    
    #Reading input data: phenotype matrix, genotype matrix, relatedness matrix 

    data = ReadData(verbose=False)

    file_pheno = resource_filename('limmbo',filename_pheno)
    file_geno = resource_filename('limmbo', filename_geno)
    file_relatedness = resource_filename('limmbo',filename_relatedness)
    
    data.getPhenotypes(file_pheno=file_pheno)
    data.getRelatedness(file_relatedness=file_relatedness,delim=",")
    data.getGenotypes(file_genotypes=file_geno)
    
    indata = InputData(verbose=False)
    indata.addPhenotypes(phenotypes = data.phenotypes)
    
    relatedness=pd.io.parsers.read_csv(file_relatedness,sep=',')
    relatedness.index = relatedness.columns
    indata.addRelatedness(relatedness = relatedness)
    
    indata.addGenotypes(genotypes=data.genotypes,genotypes_info=data.genotypes_info)
    indata.commonSamples()
    #indata.transform(transform="scale")
    
    #Getting eigenvectors and eigenvalues of relatedness matrix
    indata.S_R, indata.U_R = LA.eig(relatedness)
    
    n=indata.phenotypes.shape[0]  #sample size
    q=indata.phenotypes.shape[1]  #num phenotypes
    p=indata.genotypes.shape[1]   #num genotypes
    
    
    #computing kronecker product of identity with genotype matrix
    x= indata.genotypes.to_numpy()
    X= np.kron(np.identity(q),x)

    #vectorizing the phenotype matrix
    Y=indata.phenotypes.to_numpy()
    y=np.matrix.flatten(Y,'F')
    
   

    for i in range(maxiter):
        #print('Iteration {}'.format(i))
        #Given fixed effects, estimate covariance matrices based on standard REML
        #for larger number of phenotypes, use limmbo bootstrapped version instead
        Cg, Cn, ptime = vd_reml(indata, verbose=False)
    
        #Given covariance matrices, fit fixed effects using L1 regularization
        H=np.kron(indata.relatedness,Cg)+np.kron(np.ones((n,n)),Cn)
        H_cpu = torch.tensor(H, device='cpu')
        Hinv_cpu = torch.inverse(H_cpu)
        Hinv_cpu_s = torch.cholesky(Hinv_cpu)
        Hinv_s = Hinv_cpu_s.numpy()
    
        ytilde=Hinv_s @ y
        Xtilde=Hinv_s @ X
    
        #running regularized estination with trait-wise clustering 
        estimator = twclust(p=p, n=n, K=q)
        theta_hat = estimator.fit(Xtilde, ytilde, C1=C1, C2=C2, gamma=0.003, maxiter=1000)
        beta = estimator.Theta
    
        #update residuals for new cova estimation step
        newY=Y-x@beta
        indata.phenotypes= pd.DataFrame(newY, index=indata.pheno_samples,
                                       columns=indata.phenotype_ID)
        
    return beta, Cg, Cn
