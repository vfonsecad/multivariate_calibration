# ------------------------------------------------------------------------

#                           preprocessing functions
# by: valeria fonseca diaz
# supervisors: Wouter Saeys, Bart De Ketelaere
# ------------------------------------------------------------------------


import numpy as np
import scipy as sp
import matplotlib.pyplot as plt



# ¡¡¡ --- !!! functions to perform OSC (orthogonal signal correction)



# --------------- SIMPLS function


def simpls(xx, yy, ncp, P = None):
    
    '''
    this function simply performs the usual PLS algorithm and returns the tuple (slope, intercept, loadings)
    y = intercept + x*slope
    loadings is only in the case X scores are needed for other calculations
    
    '''
    
    assert yy.ndim == 2, 'yy needs to be a 2D array. if there is only one column, make sure it is of shape (n,1)'

    X = xx.copy()
    Y = yy.copy()
    
    N = X.shape[0]
    K = X.shape[1]
    q = Y.shape[1] 
    
    if P is None:
        P = np.identity(n = N) / N
   
    mu_x = ((P.dot(X)).sum(axis=0))/ P.sum()
    mu_y = ((P.dot(Y)).sum(axis=0))/ P.sum()


    Xc = X - mu_x
    Yc = Y - mu_y

    
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
    
    TT = Xc.dot(R)
               
        # --- final regression

        
    tcal_raw0 = np.concatenate((np.ones((X.shape[0], 1)), TT), axis=1)
    wtemp = np.linalg.solve(tcal_raw0.T.dot(P.dot(tcal_raw0)), tcal_raw0.T.dot(P.dot(Y)))          
            
    b = wtemp[0,0]    
    BPLS = R.dot(wtemp[1:, :])
 


    return (BPLS,b,R)



def osc(xx, yy, current_osc_ncp = 1, epsilon = 10e-6, total_iters = 20, osc_pls_lv = 5):

    '''
    this function performs OSC correction with the help of the simpls algorithm 
    for further information check https://www.sciencedirect.com/science/article/pii/016974399385002X
    
    output: (xx_new, (W,P,mu_x))
    
    any raw x can be then transformed as:
      
        xx_new = ((xx_raw - mu_x).dot(np.identity(n=xx_raw.shape[1]) - W.dot(np.linalg.inv(P.T.dot(W)).dot(P.T)))) + mu_x

    
    '''
    
    assert yy.ndim == 2, 'yy needs to be a 2D array. if there is only one column, make sure it is of shape (n,1)'


    X = xx.copy()
    Y = yy.copy()
    
    N = X.shape[0]
    K = X.shape[1]
    q = Y.shape[1] 

           
    A = np.identity(n = N) / N # this is here temporarily to include sample weights
    
    mu_x = ((A.dot(X)).sum(axis=0))/ A.sum()
    mu_y = ((A.dot(Y)).sum(axis=0))/ A.sum()

    Xc = X - mu_x 
    Yc = Y - mu_y

    W = np.zeros((X.shape[1],current_osc_ncp))
    P = np.zeros((X.shape[1],current_osc_ncp))
    TT = np.zeros((X.shape[0],current_osc_ncp))


    kk = 0
    while kk < current_osc_ncp:
        print('---', kk+1)

        # --- pc of xc

        xu, xs, xvt = np.linalg.svd(Xc)
        tt_old = xu[:,0:1]*xs[0]
        p = xvt.T[:,0:1]
        p = np.multiply(p, np.sign(np.sum(p)))
        
        iter_i = 0
        convergence = 10 + epsilon

        while convergence > epsilon:

            # - calculate scores
            tt = Xc.dot(p)/(p.T.dot(p))
            # - orthogonalize scores
            tt_new = (np.identity(Yc.shape[0]) - Yc.dot(np.linalg.pinv(Yc.T.dot(Yc)).dot(Yc.T))).dot(tt)
            #- update loadings
            p_new = Xc.T.dot(tt_new)/(tt_new.T.dot(tt_new))
            # - calculate convergence
            convergence = np.linalg.norm(tt_new - tt_old, axis = 0) / np.linalg.norm(tt_new, axis = 0)
            # - update scores and loadings
            tt_old = tt_new.copy()
            p = p_new.copy()

            iter_i += 1

            # - check convergence in iterations

            if iter_i > total_iters:
                convergence = 0



        # - perform regression of X and t, 5 lv by default  

        w, beta, R = simpls(Xc, tt_new, osc_pls_lv)
        w = w/np.linalg.norm(w)

        # - calculate final component that will be removed and stored

        tt = X.dot(w)
        tt = (np.identity(Yc.shape[0]) - Yc.dot(np.linalg.pinv(Yc.T.dot(Yc)).dot(Yc.T))).dot(tt)
        p = Xc.T.dot(tt)/(tt.T.dot(tt))
        Xc = Xc - tt.dot(p.T)

        # - store component

        W[:,kk] = w[:,0]
        P[:,kk] = p[:,0]
        TT[:,kk] = tt[:,0]

        kk += 1

    # --- final transformation of original data
    

    xx_new = ((xx - mu_x).dot(np.identity(n=xx.shape[1]) - W.dot(np.linalg.inv(P.T.dot(W)).dot(P.T)))) + mu_x

    return (xx_new, (W,P,mu_x))




print("finished")



# --- Perform OSC

#X_new, osc_parms = osc(X_original, Y, osc_pls_lv = 5) 

# - Note: use now X_new instead of original X in the regression model with Y. Note that OSC can be applied also after log transformations, etc. 