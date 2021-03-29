% -------------------------------------------------------------
%       create milk data prepared
%by: valeria fonseca diaz
% ------------------------------------------------------------- 
clc ; clear all;

% --- project folder and functions

mainpath = "C:/Users/u0106869/Google Drive/kul_phd/conferences/chemometrica-2019/paper/chemometrica2019_conference_paper_codes_submission02";
%mainpath = "/home/valeria/vfonsecad/kul_phd/conferences/chemometrica-2019/paper/chemometrica2019_conference_paper_codes_submission02";
addpath(mainpath + "/methods/matlab/LIBRA_20160628")
addpath(mainpath + "/methods/matlab")

% --- load data

load(mainpath + "/data/d0005_milkrobot_visnir2010/data_raw/Transmittance.mat")

y = Matrix(:,1:5);
x = Matrix(:, 10:end);

% --- robpca for split 


n = size(x,1);
all_samples = (1:1:n)';
nval = floor(n*0.2);
x_robpca = robpca(x,'k',5,'alpha',0.75);
inliers = find(x_robpca.flag.all);

KS_val = kennard_stone(x(inliers,:),nval);
val_obs = ismember(all_samples,inliers(KS_val.sample_id==1,:));


% --- get samples into cal and val

xcal = x(val_obs==0,:);
ycal = y(val_obs==0,:);
xval = x(val_obs==1,:);
yval = y(val_obs==1,:);




%% --- preprocessing Aernouts 2011

xcal_raw = xcal(:,:);
xval_raw = xval(:,:);


xcal_absorbance  = log(1./xcal_raw(:, (1000-350):(1700-350)));
xval_absorbance  = log(1./xval_raw(:, (1000-350):(1700-350)));



xcal = xcal_absorbance;
xval = xval_absorbance;





% --- save data_prepared


y_labels = char({'Fat', 'Protein', 'Lactose', 'SCC', 'Urea'});

save(mainpath + "/data/d0005_milkrobot_visnir2010/data_prepared/trans_lactose_robpca_data_prepared.mat", 'xcal','ycal', 'xval', 'yval', 'y_labels')


fprintf('------------------ FINISHED -----------------------\n')





