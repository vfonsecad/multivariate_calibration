# multivariate calibration
# by: valeria fonseca diaz valeria.fonsecadiaz@kuleuven.be
# this folder contains procedures to apply pls models for multivariate calibration

mcw-pls is an algorithm that takes simpls as a special case. this implementation is found in https://onlinelibrary.wiley.com/doi/abs/10.1002/cem.3215

## folders:

- data: data cases. usually there is a data raw folder and a data prepared folder so data is structured to be ready for analysis. a template for this is available in user documents but it is not mandatory to do it this way
- methodology: here there are the codes for mcwpls, preprocessing with osc and another procedure to read the data files as stored in the /data/d0001/dataprepared folder
- scripts: the scripts here stored make use of /data and /methodology. while the /data format is not mandatory, the methodology is. in order to use a different data input, just skip the data_class class in the first cell of the jupyter notebook, make your own data objects but keep all the rest in the first cell as it is
- user documents: there is a template to make a .mat file in matlab for the data cases




