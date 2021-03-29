 # ------------------------------------------------------------------------

#                    MCESIMPLS from paper with weights for variables
# ------------------------------------------------------------------------


from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.model_selection import cross_val_predict, KFold, GridSearchCV
import numpy as np
import matplotlib.pyplot as plt





# ------------------------- class_weighted_mcesimpls_sklearn ----------------------------------------------



class weighted_mcesimpls_sklearn(BaseEstimator, RegressorMixin):

    def __init__(self, n_components=2, max_iter=30, V_initial=None, scale_sigma2=1, robpca_check  = False, robpca_h = None):
        
        if robpca_check:
            assert max_iter == 0 and V_initial is None 
        
        self.n_components = n_components
        self.max_iter = max_iter
        self.V_initial = V_initial
        self.scale_sigma2 = scale_sigma2
        self.robpca_check = robpca_check
        self.robpca_h = robpca_h

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

                U0, sval, Qsvd = np.linalg.svd(S, full_matrices=True, compute_uv=True)
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
    
    def unimcd(self, z0, h):
        
        nz = z0.shape[0] 
        ll = nz - h + 1       
        z = np.sort(z0, axis=0)
        az = np.zeros((ll,1))
        az[0,0] = z[0:h,0].sum(axis=0)
        for samp in range(1,ll):
            az[samp,0] = az[samp-1,0] - z[samp-1,0] + z[samp+h-1] 
        az2 = np.power(az,2)/h
        sq = np.zeros((ll,1))
        sq[0,0] = np.power(z[0:h,0],2).sum(axis=0) - az2[0,0]
        for samp in range(1,ll):
            sq[samp,0] = sq[samp-1,0]-np.power(z[samp-1,0],2)+np.power(z[samp+h-1,0],2)-az2[samp,0]+az2[samp-1,0]

        sqmin = np.amin(sq, axis = 0)[0]
        ii = np.argmin(sq, axis = 0)[0]
        ndup = np.where(sq == sqmin)[0].shape[0]
        jj = np.floor((1+ndup)/2) -1
        jj = jj.astype(int)
        slutn = np.multiply(np.ones((ndup,1)),az[ii,0])
        initmean = slutn[jj,0]/h 
        initcov = sqmin / (h-1)
        # calculating consistency factor
        res = np.power(z0 - initmean,2) / initcov
        sortres = np.sort(res, axis=0)
        factor = sortres[h-1,0]/chi2.ppf(h/nz,1)
        initcov = factor*initcov
        res = np.power(z0 - initmean,2) / initcov
        quantile = chi2.ppf(0.975,1)
        weights = res<=quantile
        if not (weights.shape[0] == nz):
                weights = weights.T
        tmcd = (np.multiply(z0, weights)/weights.sum()).sum()
        smcd = np.sqrt((np.multiply(np.power(z0 - tmcd,2), weights)).sum()/(weights.sum()-1))
        
        output = {"mean" : tmcd,
                "std" : smcd}
        
        return output

        
        


    def fit(self, xx, yy):

        X = xx.copy()
        Y = yy.copy()

        N = X.shape[0]
        K = X.shape[1]
        q = Y.shape[1]
        robpca_check = self.robpca_check
        robpca_h = self.robpca_h

        P = np.identity(N)
        
        
        ncp = self.n_components
        sigma_factor = self.scale_sigma2
        max_iter = self.max_iter
        V_initial = self.V_initial
        sigma_tol = 0.000001

        if robpca_check:

            numpy2ri.activate()
            zz = np.concatenate((X, Y), axis=1)
            rospca = importr('rospca')
            current_kk = ncp + q
            my_robpca = rospca.robpca(zz, k=current_kk, h = robpca_h)
            P = np.diag(np.array(my_robpca.rx("flag.all")).flatten())
            numpy2ri.deactivate()
           
  



        
        mu_x = ((P.dot(X)).sum(axis=0)) / P.sum()
        Xc = X - mu_x
        mu_y = ((P.dot(Y)).sum(axis=0)) / P.sum()
        Yc = Y - mu_y


        if V_initial is None:

            # ---  Initial SIMPLS (RSIMPLS)

            V = self.simpls_loadings(X, Y, P, ncp)


        else:

            V = V_initial.copy()
            V.shape = (K, ncp)




        # --- Iterative process

        kk = 0
        sigma2 = 100000
        
        

        while kk < max_iter and sigma2 > sigma_tol:

            TT = Xc.dot(V)  # TT are not weighted

            tcal_raw0 = np.concatenate((np.ones((Xc.shape[0], 1)), TT), axis=1)
            wtemp = np.linalg.solve(tcal_raw0.T.dot(P.dot(tcal_raw0)), tcal_raw0.T.dot(P.dot(Y)))

            sigma_vector = np.power(Y - tcal_raw0.dot(wtemp), 2).sum(axis=1)
            sigma2 = sigma_factor * (sigma_vector).mean()
            P = np.multiply(np.exp(-sigma_vector / (2 * sigma2)), np.identity(N))

            mu_x = ((P.dot(X)).sum(axis=0)) / P.sum()
            Xc = X - mu_x
            mu_y = wtemp[0, 0]
            Yc = Y - mu_y

            V = self.simpls_loadings(X, Y, P, ncp)

            kk += 1
            
        TT = Xc.dot(V)
        P_old = np.diag(P)
        h_difference = -999698

        # --- T scores flags and final reweight for robust robpca regression
        
        if robpca_check:
            
            # --- flags for T scores
            
            
            ropca_sigmax = np.cov(X, rowvar=False, aweights=P_old)
            pp = np.linalg.solve(V.T.dot(ropca_sigmax).dot(V), V.T.dot(ropca_sigmax))    
            xtilde = TT.dot(pp);
            Cdiff = Xc-xtilde;
            P_od = np.sqrt(np.diag(Cdiff.dot(Cdiff.T)))
            P_od.shape = (N, 1)
            rank_x = np.linalg.matrix_rank(X)
            
            if not (ncp == rank_x):
                current_h = np.array(my_robpca.rx("h")).flatten()[0]
                current_h = current_h.astype(int)
                unimcd_parms = self.unimcd(np.power(P_od,2/3), h = current_h)
                cutoff_od = np.sqrt(np.power(norm.ppf(0.975, loc=unimcd_parms["mean"] , scale = unimcd_parms["std"]),3))
                
              
            else:
                
                cutoff_od = 0
                
            P_od_flag = P_od.flatten() <= cutoff_od
            
            
            # --- robpca regression
            
            zz = np.concatenate((TT, Y), axis=1)
            cov_zz = np.cov(zz, rowvar=False, aweights=P_old)
            center_zz = (np.diag(P_old).dot(zz)).sum(axis=0)/P_old.sum()
            


            rawcenterx = center_zz[0:ncp]
            rawcenterx.shape = (ncp,1)
            rawcentery = center_zz[ncp:(ncp+q)]
            rawcentery.shape = (q,1)

            rawsigmax = cov_zz[0:ncp, 0:ncp]
            rawsigmay = cov_zz[ncp:(ncp+q),ncp:(ncp+q)]

            rawsigmaxy = cov_zz[0:ncp,ncp:(ncp+q)]
            rawsigmayx = rawsigmaxy.T.copy()

            rawbeta = np.concatenate((np.linalg.solve(rawsigmax, rawsigmaxy),(rawcentery - rawsigmayx.dot(np.linalg.solve(rawsigmax, rawcenterx))).T), axis = 0)
            
            
            rawcovE = rawsigmay - (rawbeta[0:ncp,0:q].T.dot(rawsigmax).dot(rawbeta[0:ncp,0:q]))
            X_extended = np.concatenate((TT,np.ones((TT.shape[0], 1))), axis=1)

            rawfitted = X_extended.dot(rawbeta)
            residuals = Y - rawfitted
            rawcovE_inv = np.linalg.inv(rawcovE)

            
            P_resd_flag = np.sqrt(np.diag(residuals.dot(rawcovE_inv).dot(residuals.T))) <= np.sqrt(chi2.ppf(0.975,q))


            P_old_2 = P_resd_flag*1

            # -- second stange robpca regression

            cov_zz = np.cov(zz, rowvar=False, aweights=P_old_2)
            center_zz = (np.diag(P_old_2).dot(zz)).sum(axis=0)/P_old_2.sum()



            rawcenterx = center_zz[0:ncp]
            rawcenterx.shape = (ncp,1)
            rawcentery = center_zz[ncp:(ncp+q)]
            rawcentery.shape = (q,1)

            rawsigmax = cov_zz[0:ncp, 0:ncp]
            rawsigmay = cov_zz[ncp:(ncp+q),ncp:(ncp+q)]

            rawsigmaxy = cov_zz[0:ncp,ncp:(ncp+q)]
            rawsigmayx = rawsigmaxy.T.copy()

            rawbeta = np.concatenate((np.linalg.solve(rawsigmax, rawsigmaxy),(rawcentery - rawsigmayx.dot(np.linalg.solve(rawsigmax, rawcenterx))).T), axis = 0)


            rawcovE = rawsigmay - (rawbeta[0:ncp,0:q].T.dot(rawsigmax).dot(rawbeta[0:ncp,0:q]))
            X_extended = np.concatenate((TT,np.ones((TT.shape[0], 1))), axis=1)

            rawfitted = X_extended.dot(rawbeta)
            residuals = Y - rawfitted
            rawcovE_inv = np.linalg.inv(rawcovE)


            P_resd_flag_2 = np.sqrt(np.diag(residuals.dot(rawcovE_inv).dot(residuals.T))) <= np.sqrt(chi2.ppf(0.975,q))



            P = np.diag(1*np.multiply(P_od_flag, P_resd_flag))
            
                    
            
            
            
            
        # --- wmcesimpls final regression

        
        tcal_raw0 = np.concatenate((np.ones((X.shape[0], 1)), TT), axis=1)
        wtemp = np.linalg.solve(tcal_raw0.T.dot(P.dot(tcal_raw0)), tcal_raw0.T.dot(P.dot(Y)))
            
            
        
        BPLS = V.dot(wtemp[1:, :])

        self.x_scores_ = TT
        self.x_weights_ = V
        self.x_mu = mu_x
        self.y_mu = wtemp[0, 0]
        self.BPLS = BPLS
        self.sample_weights = P
        self.x_scores_coef_ = wtemp[1:,:]

    def predict(self, X):

        Ypred = self.y_mu + (X - self.x_mu).dot(self.BPLS)

        return Ypred


# ------------------------- class_weighted_mcesimpls ----------------------------------------------


class weighted_mcesimpls(object):

    def __init__(self, xx, yy, n_components, model_name=""):
        ''' Initialize a Weighted mcesimpls proposal Class object with calibration data '''

        assert type(xx) is np.ndarray and type(yy) is np.ndarray and xx.shape[0] == yy.shape[0]

        self.xcal = xx.copy()
        self.ycal = yy.copy()
        self.Ncal = xx.shape[0]
        self.XK = xx.shape[1]
        self.YK = yy.shape[1]
        self.model_name = model_name
        self.n_components = n_components

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

    def train(self, iters=30, current_v0=None, factor_sigma=1, robpca_check = False, robpca_h = 30):

        ''' robpca will make use of the rospca library in R. This is mean to be True when RSIMPLS is to be performed'''

        weighted_mcesimpls_trainObject = weighted_mcesimpls_sklearn(n_components=self.n_components, max_iter=iters, V_initial=current_v0,
                                             scale_sigma2=factor_sigma, robpca_check = robpca_check, robpca_h = robpca_h)
        weighted_mcesimpls_trainObject.fit(self.get_xcal(), self.get_ycal())

        weighted_mcesimpls_fitted = weighted_mcesimpls_trainObject.predict(self.get_xcal())
        weighted_mcesimpls_coeff = weighted_mcesimpls_trainObject.BPLS


        weighted_mcesimpls_coeff.shape = (self.XK, self.YK)
        weighted_mcesimpls_fitted.shape = (self.Ncal, self.YK)

        weighted_mcesimpls_Output = {'BPLS': weighted_mcesimpls_coeff,
                         'x_mean':weighted_mcesimpls_trainObject.x_mu,
                         'y_mean':weighted_mcesimpls_trainObject.y_mu,
                         'x_scores':weighted_mcesimpls_trainObject.x_scores_,
                         'x_weights':weighted_mcesimpls_trainObject.x_weights_,
                         'x_scores_coef' : weighted_mcesimpls_trainObject.x_scores_coef_,
                         'fitted': weighted_mcesimpls_fitted,
                         'trainObject': weighted_mcesimpls_trainObject,
                         'sample_weights': weighted_mcesimpls_trainObject.sample_weights,
                         'factor_sigma': factor_sigma,
                         'rsimpls' : robpca_check
                         

                         }

        return weighted_mcesimpls_Output


        # --- CrossValidation

    def crossval_KFold(self, trainObject, number_splits=10):

        cvObject = KFold(n_splits=number_splits)
        cv_predicted = cross_val_predict(trainObject, self.get_xcal(), self.get_ycal(), cv=cvObject)

        cv_Output = {'cvPredicted': cv_predicted}

        return cv_Output


    def tune_sigma_factor(self, sigma_factor_range,robpca_check=False):



        weighted_mcesimpls_trainObject = weighted_mcesimpls_sklearn(n_components=self.n_components, max_iter=30, V_initial=None, robpca_check=robpca_check)
        cv_sigma_factor = {'scale_sigma2': list(sigma_factor_range)}
        cvObject = KFold(n_splits=10)

        TuneCV = GridSearchCV(estimator=weighted_mcesimpls_trainObject, param_grid=cv_sigma_factor, cv=cvObject,
                              scoring='neg_mean_squared_error', return_train_score=True)
        TuneCV.fit(self.get_xcal(), self.get_ycal())

        tune_Output = {'rmsecv': np.sqrt(-1 * TuneCV.cv_results_['mean_train_score']),
                       'grid': sigma_factor_range}

        return tune_Output

    def predict(self, X, weighted_mcesimpls_Output):

        Ypred = weighted_mcesimpls_Output["y_mean"] + (X - weighted_mcesimpls_Output["x_mean"]).dot(weighted_mcesimpls_Output["BPLS"])

        return Ypred

    def predict_x(self, X, weighted_mcesimpls_Output):

        Xpred = (X - weighted_mcesimpls_Output["x_mean"]).dot(weighted_mcesimpls_Output["x_weights"])

        return Xpred
