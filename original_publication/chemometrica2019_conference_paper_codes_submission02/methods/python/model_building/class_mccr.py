 # ------------------------------------------------------------------------

#                    MCCR from paper 
# ------------------------------------------------------------------------


from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.model_selection import cross_val_predict, KFold, GridSearchCV
import numpy as np
import matplotlib.pyplot as plt





# ------------------------- class_weighted_mcesimpls_sklearn ----------------------------------------------



class mccr_sklearn(BaseEstimator, RegressorMixin):

    def __init__(self, reg_lambda=2, max_iter=30, scale_sigma2=1):
        
        
        self.reg_lambda = reg_lambda
        self.max_iter = max_iter
        self.scale_sigma2 = scale_sigma2

  

    def fit(self, xx, yy):

        X = xx.copy()
        Y = yy.copy()

        N = X.shape[0]
        K = X.shape[1]
        q = Y.shape[1]
      
        P = np.identity(N)
        
        
        reg_lambda = self.reg_lambda
        sigma_factor = self.scale_sigma2
        max_iter = self.max_iter
      
        sigma_tol = 0.000001    

        
        mu_x = ((P.dot(X)).sum(axis=0)) / P.sum()
        Xc = X - mu_x
             
        
       
        Xreg = np.concatenate((np.ones((N,1)),Xc), axis = 1)
        V = np.linalg.solve(Xreg.T.dot(Xreg) + reg_lambda*np.eye(K+1), Xreg.T.dot(Y))
        V.shape = (K+1, q)
        

        # --- Iterative process

        kk = 0
        sigma2 = 100000
        
        

        while kk < max_iter and sigma2 > sigma_tol:           

            sigma_vector = np.power(Y - Xreg.dot(V), 2).sum(axis=1)
            sigma2 = sigma_factor * (sigma_vector).mean()
            P = np.multiply(np.exp(-sigma_vector / (2 * sigma2)), np.identity(N))

            mu_x = ((P.dot(X)).sum(axis=0)) / P.sum()
            Xc = X - mu_x
            
            
            Xreg = np.concatenate((np.ones((N,1)),Xc), axis = 1)
            V = np.linalg.solve(Xreg.T.dot(P).dot(Xreg) + reg_lambda*np.eye(K+1), Xreg.T.dot(P).dot(Y))
            V.shape = (K+1, q)
           
            kk += 1
            
            
        Vtemp = V[1:, :]
        TT = Xc.dot(Vtemp)
       
            
        # --- mccr final regression

        
        tcal_raw0 = np.concatenate((np.ones((X.shape[0], 1)), TT), axis=1)
        wtemp = np.linalg.solve(tcal_raw0.T.dot(P.dot(tcal_raw0)), tcal_raw0.T.dot(P.dot(Y)))
            
            
        
        BPLS = Vtemp.dot(wtemp[1:, :])

        self.x_scores_ = TT
        self.x_weights_ = Vtemp
        self.x_mu = mu_x
        self.y_mu = wtemp[0, 0]
        self.BPLS = BPLS
        self.sample_weights = P
        self.x_scores_coef_ = wtemp[1:,:]

    def predict(self, X):

        Ypred = self.y_mu + (X - self.x_mu).dot(self.BPLS)

        return Ypred


# ------------------------- class_weighted_mcesimpls ----------------------------------------------


class mccr(object):

    def __init__(self, xx, yy, reg_lambda, model_name=""):
        ''' Initialize a MCCR class object with calibration data '''

        assert type(xx) is np.ndarray and type(yy) is np.ndarray and xx.shape[0] == yy.shape[0]

        self.xcal = xx.copy()
        self.ycal = yy.copy()
        self.Ncal = xx.shape[0]
        self.XK = xx.shape[1]
        self.YK = yy.shape[1]
        self.model_name = model_name
        self.reg_lambda = reg_lambda

    def __str__(self):
        return 'class_mccr'

        # --- Copy of data

    def get_xcal(self):
        ''' Get copy of xcal data '''
        return self.xcal

    def get_ycal(self):
        ''' Get copy of ycal data '''
        return self.ycal

       # --- Define performance measures

    def rmse(self, yy, y_pred, sample_weights=None):

        if sample_weights is None:
            sample_weights = np.ones((yy.shape[0], 1)) / yy.shape[0]
        else:
            sample_weights = sample_weights / sample_weights.sum(axis=0)

        r = yy.copy() - y_pred.copy()
        r = r ** 2
        msep = np.average(r, axis=0, weights=sample_weights)
        rmse = np.sqrt(msep)

        return rmse
    
    def r2(self, yy, y_pred, sample_weights=None):

        if sample_weights is None:
            sample_weights = np.ones((yy.shape[0], 1)) / yy.shape[0]
        else:
            sample_weights = sample_weights / sample_weights.sum(axis=0)
            
        
        P = np.diag(sample_weights.flatten())
        yy_reg = np.concatenate((np.ones((yy.shape[0], 1)), yy), axis=1)
        coeffs = np.linalg.solve(yy_reg.T.dot(P.dot(yy_reg)), yy_reg.T.dot(P.dot(y_pred)))
        y_pred_fitted = yy_reg.dot(coeffs)
        
        mu_y_pred = ((P.dot(y_pred)).sum(axis=0)) / P.sum()
        y_pred_c = y_pred - mu_y_pred
        y_pred_res = y_pred_fitted - y_pred
        total_ss = ((P.dot(y_pred_c**2)).sum(axis=0)) 
        residual_ss = ((P.dot(y_pred_res**2)).sum(axis=0))
        
        r2 = 1 - (residual_ss / total_ss)


        return r2

        # --- Define Weighted mcesimpls regression from my class in scikit learn

    def train(self, iters=30, factor_sigma=1):

        ''' '''

        mccr_trainObject = mccr_sklearn(reg_lambda=self.reg_lambda, max_iter=iters, scale_sigma2=factor_sigma)
        mccr_trainObject.fit(self.get_xcal(), self.get_ycal())

        mccr_fitted = mccr_trainObject.predict(self.get_xcal())
        mccr_coeff = mccr_trainObject.BPLS


        mccr_coeff.shape = (self.XK, self.YK)
        mccr_fitted.shape = (self.Ncal, self.YK)

        mccr_Output = {'BPLS': mccr_coeff,
                         'x_mean':mccr_trainObject.x_mu,
                         'y_mean':mccr_trainObject.y_mu,
                         'x_scores':mccr_trainObject.x_scores_,
                         'x_weights':mccr_trainObject.x_weights_,
                         'x_scores_coef' : mccr_trainObject.x_scores_coef_,
                         'fitted': mccr_fitted,
                         'trainObject': mccr_trainObject,
                         'sample_weights': mccr_trainObject.sample_weights,
                         'factor_sigma': factor_sigma,
                         'lambda' : self.reg_lambda
                         

                         }

        return mccr_Output


        # --- CrossValidation

    def crossval_KFold(self, trainObject, number_splits=10):

        cvObject = KFold(n_splits=number_splits)
        cv_predicted = cross_val_predict(trainObject, self.get_xcal(), self.get_ycal(), cv=cvObject)

        cv_Output = {'cvPredicted': cv_predicted}

        return cv_Output




    def predict(self, X, mccr_Output):

        Ypred = mccr_Output["y_mean"] + (X - mccr_Output["x_mean"]).dot(mccr_Output["BPLS"])

        return Ypred

