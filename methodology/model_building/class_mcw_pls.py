# ------------------------------------------------------------------------

#                               mcw-pls
# by: valeria fonseca diaz
# supervisors: Wouter Saeys, Bart De Ketelaere
# ------------------------------------------------------------------------


from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.model_selection import cross_val_predict, KFold, GridSearchCV
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt



# ------------------------- mcw_pls_sklearn ----------------------------------------------



class mcw_pls_sklearn(BaseEstimator, RegressorMixin):

    def __init__(self, n_components=2, max_iter=30, R_initial=None, scale_sigma2=1):

        
        self.n_components = n_components
        self.max_iter = max_iter
        self.R_initial = R_initial
        self.scale_sigma2 = scale_sigma2
    
    def simpls_loadings(self,xx, yy, P, ncp):
        
        
        
        mu_x = ((P.dot(xx)).sum(axis=0))/ P.sum()
        mu_y = ((P.dot(yy)).sum(axis=0))/ P.sum()
        

        Xc = xx.copy() - mu_x
        Yc = yy.copy() - mu_y

        N = Xc.shape[0]
        K = Xc.shape[1]
        q = Yc.shape[1]

        R = np.zeros((K, ncp))  # Weights to get T components
        V = np.zeros((K, ncp))  # orthonormal base for X loadings
        S = Xc.T.dot(P).dot(Yc)  # cov matrix

        aa = 0

        while aa < ncp:
            
            r = S[:,:]            
            
            if q > 1:

                U0, sval, Qsvd = sp.linalg.svd(S, full_matrices=True, compute_uv=True)
                Sval = np.zeros((U0.shape[0], Qsvd.shape[0]))
                Sval[0:sval.shape[0], 0:sval.shape[0]] = np.diag(sval)

                Rsvd = U0.dot(Sval)
                r = Rsvd[:, 0]
                
                
            tt = Xc.dot(r)
            tt.shape = (N, 1)
            tt = tt - ((P.dot(tt)).sum(axis=0)/ P.sum())
            TT_scale = np.sqrt(tt.T.dot(P).dot(tt))
            # Normalize
            tt = tt / TT_scale            
            r = r / TT_scale
            r.shape = (K, 1)
            p = Xc.T.dot(P).dot(tt)
            v = p
            v.shape = (K, 1)

            if aa > 0:
                v = v - V.dot(V.T.dot(p))

            v = v / np.sqrt(v.T.dot(v))
            S = S - v.dot(v.T.dot(S))

            R[:, aa] = r[:, 0]
            V[:, aa] = v[:, 0]

            aa += 1

        return R
    


    
    def fit(self, xx, yy):

        X = xx.copy()
        Y = yy.copy()

        N = X.shape[0]
        K = X.shape[1]
        q = Y.shape[1]
        
      
            
        P = np.identity(N)
        
        
        ncp = self.n_components
        sigma_factor = self.scale_sigma2
        max_iter = self.max_iter
        R_initial = self.R_initial
        sigma_tol = 0.000001

  

        
        mu_x = ((P.dot(X)).sum(axis=0)) / P.sum()
        Xc = X - mu_x
        mu_y = ((P.dot(Y)).sum(axis=0)) / P.sum()
        Yc = Y - mu_y


        if R_initial is None:

            # ---  Initial SIMPLS (RSIMPLS)

            R = self.simpls_loadings(X, Y, P, ncp)


        else:

            R = R_initial.copy()
            R.shape = (K, ncp)




        # --- Iterative process

        kk = 0
        sigma2 = 100000
        
        

        while kk < max_iter and sigma2 > sigma_tol:

            TT = Xc.dot(R)  # TT are not weighted

            tcal_raw0 = np.concatenate((np.ones((Xc.shape[0], 1)), TT), axis=1)
            wtemp = sp.linalg.solve(tcal_raw0.T.dot(P.dot(tcal_raw0)), tcal_raw0.T.dot(P.dot(Y)))

            sigma_vector = np.power(Y - tcal_raw0.dot(wtemp), 2).sum(axis=1)
            sigma2 = sigma_factor * (sigma_vector).mean()
            P = np.multiply(np.exp(-sigma_vector / (2 * sigma2)), np.identity(N))

            mu_x = ((P.dot(X)).sum(axis=0)) / P.sum()
            Xc = X - mu_x
            mu_y = wtemp[0:1, :]
            Yc = Y - mu_y

            R = self.simpls_loadings(X, Y, P, ncp)

            kk += 1
            
        TT = Xc.dot(R)
               
        # --- mcw-pls final regression

        
        tcal_raw0 = np.concatenate((np.ones((X.shape[0], 1)), TT), axis=1)
        wtemp = sp.linalg.solve(tcal_raw0.T.dot(P.dot(tcal_raw0)), tcal_raw0.T.dot(P.dot(Y)))
            
            
        
        BPLS = R.dot(wtemp[1:, :])

        self.x_scores_ = TT
        self.x_weights_ = R
        self.x_mu = mu_x
        self.y_mu = wtemp[0:1, :]
        self.BPLS = BPLS
        self.sample_weights = P
        self.x_scores_coef_ = wtemp[1:,:]

    def predict(self, X):

        Ypred = self.y_mu + (X - self.x_mu).dot(self.BPLS)

        return Ypred


# ------------------------- mcw_pls ----------------------------------------------


class mcw_pls(object):

    def __init__(self, xx, yy, lv, model_name=""):
        ''' Initialize a Weighted mcesimpls proposal Class object with calibration data '''

        assert type(xx) is np.ndarray and type(yy) is np.ndarray and xx.shape[0] == yy.shape[0]

        self.xcal = xx.copy()
        self.ycal = yy.copy()
        self.Ncal = xx.shape[0]
        self.XK = xx.shape[1]
        self.YK = yy.shape[1]
        self.model_name = model_name
        self.lv = lv

    def __str__(self):
        return 'class_weighted_mcesimpls'

        # --- Copy of data

    def get_xcal(self):
        ''' Get copy of xcal data '''
        return self.xcal

    def get_ycal(self):
        ''' Get copy of ycal data '''
        return self.ycal

        # --- Define y names and wv numbers (in nm)

    def set_wv_varlabel(self, x_wv0):
        x_wv = np.array(x_wv0).flatten()
        self.wv_varlabel = x_wv

    def set_yy_varlabel(self, y_constituent0):
        y_names = np.array(y_constituent0).flatten()
        self.yy_varlabel = y_names

    def plot_spectra(self, xx):
        fig, ax = plt.subplots()
        ax.plot(xx.T)
        ticks = np.round(np.arange(0, self.XK, self.XK / 6))
        plt.xticks(ticks, np.round(self.wv_varlabel[ticks.astype(int)].flatten(), 1))
        plt.xlabel("Wavelength (nm)")
        plt.title(self.model_name)
        plt.show()

       # --- Define performance measures

    def rmse(self, yy, y_pred, sample_weights=None):
        
        
        N = yy.shape[0]

        if sample_weights is None:
            P = np.identity(n = N) / N
        else:
            sample_weights_vec = sample_weights.flatten()
            P = np.diag(v = sample_weights_vec) / sample_weights_vec.sum()

    
        r = yy.copy() - y_pred.copy()
        r = np.power(r,2)
        msep = ((P.dot(r)).sum(axis=0)) / P.sum()
        rmse = np.sqrt(msep)

        return rmse
    


        # --- define mcw-pls regression from my class in scikit learn

    def train(self, iters=30, current_R0=None, factor_sigma=1):
        
        '''
        'BPLS': regression vector
        'x_mean': mean of x
                         'y_mean':mean of y
                         'x_scores': T scores
                         'x_weights': x loadings
                         'x_scores_coef' : coefficients of regression between T scores and Y
                         'fitted': fitted y values
                         'train_object': sklearn train object,
                         'sample_weights': sample weights,
                         'factor_sigma': sigma factor
        
        '''

        

        mcw_pls_train_object = mcw_pls_sklearn(n_components=self.lv, max_iter=iters, R_initial=current_R0,scale_sigma2=factor_sigma)
        mcw_pls_train_object.fit(self.get_xcal(), self.get_ycal())

        mcw_pls_fitted = mcw_pls_train_object.predict(self.get_xcal())
        mcw_pls_coeff = mcw_pls_train_object.BPLS


        mcw_pls_coeff.shape = (self.XK, self.YK)
        mcw_pls_fitted.shape = (self.Ncal, self.YK)

        mcw_pls_output = {'BPLS': mcw_pls_coeff,
                         'x_mean':mcw_pls_train_object.x_mu,
                         'y_mean':mcw_pls_train_object.y_mu,
                         'x_scores':mcw_pls_train_object.x_scores_,
                         'x_weights':mcw_pls_train_object.x_weights_,
                         'x_scores_coef' : mcw_pls_train_object.x_scores_coef_,
                         'fitted': mcw_pls_fitted,
                         'train_object': mcw_pls_train_object,
                         'sample_weights': mcw_pls_train_object.sample_weights,
                         'factor_sigma': factor_sigma
                         }

        return mcw_pls_output


        # --- cross-validation

    def crossval_KFold(self, train_object, number_splits=10):

        cv_object = KFold(n_splits=number_splits)
        cv_predicted = cross_val_predict(train_object, self.get_xcal(), self.get_ycal(), cv=cv_object)

        cv_output = {'cv_predicted': cv_predicted}

        return cv_output



    def predict(self, X, mcw_pls_output):

        y_pred = mcw_pls_output["y_mean"] + (X - mcw_pls_output["x_mean"]).dot(mcw_pls_output["BPLS"])

        return y_pred

    def predict_x(self, X, mcw_pls_output):

        x_pred = (X - mcw_pls_output["x_mean"]).dot(mcw_pls_output["x_weights"])

        return x_pred


