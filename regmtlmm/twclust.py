#code for trait-wise clustering

import numpy as np
from sklearn.metrics import mean_squared_error
#import glmnet_python
#from glmnet import glmnet, glmnetSet

from sklearn import linear_model

class twclust(object):
    
    def __init__(self, p=20, n=80, K=15):
        
        self.p, self.n, self.K = p, n, K
        
        
    def fit(self, X, Y, C1=1000, C2=20, gamma=0.003, maxiter=1000):
        theta_hat = self.individual_lasso(X, Y) #compute lasso estimator
        WW, L = self.give_weights(theta_hat) #based on estimator determine weights indicating similarity between phenotype pairs
        
        Theta = self.prox_decomp(X, Y, WW, L, C1=C1, C2=C2, gamma=gamma, maxiter=maxiter)
        
        self.Theta = Theta
        self.WW = WW
        self.L = L
        
        return self
    
    
    def predict(self, X_test):
        
        return X_test @ self.Theta
    
    
    def RMSE(self, Y_true, Y_pred):
        
        return (mean_squared_error(Y_true, Y_pred))**0.5
    
    #compute lasso estimation  
    def individual_lasso(self, X, Y):
        
        clf=linear_model.Lasso(alpha=0.2,fit_intercept=False)
        clf.fit(X=X,y=Y)
    
        Theta_hat=np.reshape(clf.coef_,(self.p,self.K),'F')
     
        
        return Theta_hat
    
    #Based lasso estimation, compute weights to determine which pairs of traits should be encouraged to be clustered together 
    def give_weights(self, Theta_hat, k=5):
        
        YY = Theta_hat.copy()
        YY = YY-np.mean(YY)
        YY = YY/np.linalg.norm(YY)
        
        D = np.ones([self.K, self.K]) * 1e10
        for i in range(self.K):
            for j in range(self.K):
                if i!=j:
                    D[i, j] = np.linalg.norm(YY[:, i]-YY[:, j])
                    
        Near = np.zeros([self.K, k])
        for i in range(self.K):
            d = D[i, :]
            I = np.argsort(d)
            Near[i, :] = I[:k]
        Near = np.array(Near, dtype=int)
            
        phi = 30
        W = np.zeros([self.K, self.K])
        for i in range(self.K):
            for j in Near[i]:
                if i < j:
                    W[i,j] = np.exp( -phi * np.linalg.norm(YY[:,i]-YY[:,j])**2 )
                else:
                    W[j,i] = np.exp( -phi * np.linalg.norm(YY[:,i]-YY[:,j])**2 )
                    
        WW = W + W.T
        
        L = []
        for i in range(self.K):
            for j in Near[i]:
                L.append(np.sort([i, j]))
        L = np.unique(np.array(L), axis=0)
                
        return WW, L
    
    #perform estimation with trait-wise clustering
    def prox_decomp(self, X, Y, WW, L, C1=10, C2=10, gamma=0.003, maxiter=1000):
        
        p, n, K = self.p, self.n, self.K
        
        lambda1 = np.sqrt(np.log(p)/n) * C1
        lambda2 = np.sqrt(np.log(p)/n) * C2
        
        m = L.shape[0] + 2
        w = np.ones([m,1])/m
        y = np.random.randn(p*K,m)
        x0 = y @ w
        x = x0
        
        iteration = 0
        distance = 1
        
        while distance > 1e-5:
            if iteration > maxiter:
                break
                
            iteration += 1
            x_old = x
            P = np.zeros([p*K,m])
            
            # f1
            y1 = y[:,0]
            sig = gamma/w[0]
            p1 = np.linalg.inv(sig*X.T@X + 0.5*np.identity(p*K)) @ (sig*X.T@Y + 0.5 * y1)
            P[:, 0] = p1
            
            # sparsity
            y2 = y[:,1]
            sig = gamma/w[1]
            s = 1 - (lambda1*sig) / np.absolute(y2)
            P[:, 1] = np.maximum(s,0) * y2
            
            # clustering
            for ii in range(2, m):
                pair = L[ii-2,:]
                yii = y[:,ii]
                i = pair[0] 
                j = pair[1]
                # i_start = i*p
                # i_end = (i+1)*p
                sig = gamma/w[ii]
                P[:,ii] = yii
                con = sig * lambda2 * WW[i,j] / np.linalg.norm(yii[i*p:(i+1)*p] - yii[j*p:(j+1)*p])
                P[i*p:(i+1)*p,ii] = (1-con) * yii[i*p:(i+1)*p] + con * yii[j*p:(j+1)*p]
                P[j*p:(j+1)*p,ii] = con * yii[i*p:(i+1)*p] + (1-con) * yii[j*p:(j+1)*p]
            
            
            
            Pn = P @ w
            lambdan = 1

            y = y + lambdan * ( (2*Pn - x) @ np.ones([1,m]) - P )
            x = x + lambdan * (Pn - x)
    
            distance = np.linalg.norm(x - x_old)
            
            if iteration % 1 == 0:
                # print('Iteration {}: {}'.format(iteration, distance))
                pass
            
        Theta_splitting = np.reshape(x, [p,K])
        
        #print('Iteration {}: {}'.format(iteration, distance))
        
        return Theta_splitting
            
            
        
        
        
        
        
        
            
        
            
        
        