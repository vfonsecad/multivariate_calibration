# ----------------------------------------------------------------------

# --- Class import data for chemometrics analysis

# ----------------------------------------------------------------------


import scipy.io as sp_io
import numpy as np

class chemometrics_data(object):
    
    def __init__(self,  mat_filename, 
                        include_cal = True, 
                        include_val = True,
                        include_test = False, 
                        include_unlabeled = False,
                        y_all_range = True , x_all_range = True,
                        obs_all_cal = True, obs_all_val = True,
                        obs_all_test = True, obs_all_unlabeled = True,
                        y_range = None, x_range = None,
                        obs_cal = None, obs_val = None, obs_test = None,
                        obs_unlabeled = None,
                        data_identifier = None,
                        shuffle = False):
        
        
    
        '''
        Read chemometrics data for calibration models
        
        
        mat_filename: name of the matfile that contains all elements for the chemometrics data 
        include_cal: Default True. 'xcal' and 'ycal' exist in the matfile
        include_val: Default True. 'xval' and 'yval' exist in the matfile
        include_test: Default True. 'xtest' and 'ytest' exist in the matfile
        include_unlabeled: Default False. 'x_unlabeled' and 'y_unlabeled' exist in the matfile
        y_all_range : Default True. Include all columns (reference values) of Y
        x_all_range : Default True. Include all columns (wavelengths) of X
        obs_all_cal: Default True. Include all cal samples. 
        obs_all_val : Default True. Include all val samples. Not used if include_val is False 
        obs_all_test : Default True. Include all test samples. Not used if include_test is False
        obs_all_unlabeled : Default True. Include all unlabeled  samples. Not used if include_unlabeled is False
        y_range : Default None. 1-D numpy array of Y columns to select. If None, all columns are selected
        x_range: Default None. 1-D numpy array of X columns to select. If None, all columns are selected
        obs_cal: Default None. 1-D numpy array of xcal rows to select. If None, all rows are selected
        obs_val : Default None. 1-D numpy array of xval rows to select. If None, all rows are selected
        obs_test : Default None. 1-D numpy array of xtest rows to select. If None, all rows are selected
        obs_unlabeled : Default None. 1-D numpy array of x_unlabeled rows to select. If None, all rows are selected
        data_identifier : string of a name to identify this data case

        '''
        
        assert include_cal is True 
        assert data_identifier is not None and type(data_identifier) is str
        
        if not y_all_range:
            assert type(y_range) is np.ndarray 
        if not x_all_range:
            assert type(x_range) is np.ndarray
        if not obs_all_cal:
            assert type(obs_cal) is np.ndarray
        if not obs_all_val:
            assert type(obs_val) is np.ndarray
        if not obs_all_test:
            assert type(obs_test) is np.ndarray
        if not obs_all_unlabeled:
            assert type(obs_unlabeled) is np.ndarray
            
            


        # --- Read Data

        data_mat = sp_io.loadmat(mat_filename)


            # --- Calibration data and basic settings

        ycal_00 = data_mat['ycal'][:,:]
        xcal_00 = data_mat['xcal'][:,:]
        ncal = xcal_00.shape[0]
        K = xcal_00.shape[1]
        YK = ycal_00.shape[1]
        
        # --- Select columns
        
        if y_all_range:            
            y_range = np.arange(0, YK)
        if x_all_range:
            x_range = np.arange(0, K)
        
        ycal_01 = ycal_00[:, y_range]
        xcal_01 = xcal_00[:, x_range]
        
        # --- Select rows
        
        if obs_all_cal:             
            obs_cal = np.arange(0, ncal)
        
        if shuffle:
            obs_cal = np.random.permutation(obs_cal)
        
        ycal_02 = ycal_01[obs_cal, :]
        xcal_02 = xcal_01[obs_cal, :]
        
        y_names = data_mat["y_labels"][y_range]
        data_included = []
        data_included.append("cal")
        ncal = ycal_02.shape[0]
        
        self.ncal = ncal
        self.K = K
        self.YK = YK
        self.ycal = ycal_02
        self.xcal = xcal_02
        self.y_names = y_names
        self.data_identifier = data_identifier + "!*" + "-".join(list(y_names)) + "*!"
        
       
        
            # --- Validation data

        if include_val: 
            
            data_included.append("val")

            yval_01 = data_mat['yval'][:,y_range]
            xval_01 = data_mat['xval'][:,x_range]
            nval = xval_01.shape[0]
            
                # --- Select rows

            if obs_all_val:             
                obs_val = np.arange(0, nval)

            yval_02 = yval_01[obs_val, :]
            xval_02 = xval_01[obs_val, :]
            
            self.nval = nval
            self.yval = yval_02
            self.xval = xval_02

            # --- Test data

        if include_test:   
            
            data_included.append("test")

            ytest_01 = data_mat['ytest'][:,y_range]
            xtest_01 = data_mat['xtest'][:,x_range]
            ntest = xtest_01.shape[0]
        
            # --- Select rows

            if obs_all_test:             
                obs_test = np.arange(0, ntest)

            ytest_02 = ytest_01[obs_test, :]
            xtest_02 = xtest_01[obs_test, :]
            
            self.ntest = ntest
            self.ytest = ytest_02
            self.xtest = xtest_02

            # --- Unlabeled data

        if include_unlabeled:  
            
            data_included.append("unlabeled")
           
            x_unlabeled_01 = data_mat['x_unlabeled'][:,x_range]
            n_unlabeled = x_unlabeled_01.shape[0]
            
             # --- Select rows

            if obs_all_unlabeled:             
                obs_unlabeled = np.arange(0, n_unlabeled)

            x_unlabeled_02 = x_unlabeled_01[obs_unlabeled, :]
            
            self.n_unlabeled = n_unlabeled
            self.x_unlabeled = x_unlabeled_02
            
        self.data_included = data_included
        
    def get_cal(self):
        
        return {"ycal": self.ycal,
                "xcal": self.xcal}
    
    def get_val(self):
        
        return {"yval": self.yval,
                "xval": self.xval}
    
    def get_test(self):
        
        return {"ytest": self.ytest,
                "xtest": self.xtest}
        
    def get_unlabeled(self):
        
        return {"x_unlabeled": self.x_unlabeled}
        
            

