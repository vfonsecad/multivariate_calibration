% -------------------------------------------------------------
%       create corn data prepared
%by: valeria fonseca diaz
% ------------------------------------------------------------- 
clc ; clear all;

% --- project folder and functions

mainpath = "/home/valeria/vfonsecad/kul_phd/conferences/chemometrica-2019/paper/chemometrica2019_conference_paper_codes_submission02";
addpath(mainpath + "/methods/matlab/LIBRA_20160628")
addpath(mainpath + "/methods/matlab")

% --- load data

load(mainpath + "/data/d0001_corn/data_raw/corn.mat")
x = mp6spec.data;
y = propvals.data;

% --- robpca and kennard stone to select val samples

rng(0)
x_robpca = robpca(x,'kmax',15,'alpha',0.75);
all_samples = (1:size(x,1))';
inliers = find(x_robpca.flag.all);

KS_val = kennard_stone(x(inliers,:),16);
val_obs = ismember(all_samples,inliers(KS_val.sample_id==1,:));


% --- get samples into cal and val

xcal = x(val_obs==0,:);
ycal = y(val_obs==0,:);
xval = x(val_obs==1,:);
yval = y(val_obs==1,:);



% --- save data_prepared

y_labels = propvals.label{2};

save(mainpath + "/data/d0001_corn/data_prepared/corn_data_prepared.mat", 'xcal','ycal', 'xval', 'yval', 'y_labels')





