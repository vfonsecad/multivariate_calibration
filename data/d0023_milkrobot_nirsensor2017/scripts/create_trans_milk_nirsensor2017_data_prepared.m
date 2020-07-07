% -------------------------------------------------------------
%       create milk d0023 data prepared
%by: valeria fonseca diaz
% ------------------------------------------------------------- 
clc ; clear all;

% --- project folder and functions

mainpath = "/home/valeria/vfonsecad/kul_phd/programming/phd_valeria_fonseca_diaz_wp1/wp1_study001_robust_multicalibration";
addpath(mainpath + "/methodology/matlab/level1/LIBRA_20160628_NEWVERSION")


% --- load data

load(mainpath + "/data/d0023_milkrobot_nirsensor2017/data_raw/Data_Tot.mat")

xdata = Trans_Tot;
ydata = [Cow_ID_Tot, SET, Fat_Tot, Lact_Tot, Prot_Tot, Urea_Tot];

keep_obs = sum(ismissing(ydata),2) == 0;

xdata_complete = xdata(keep_obs,:);
ydata_complete = ydata(keep_obs,:);

% dataset 1 (weeks 1 to 5 for calibration, 6 for validation and 7-8 for test)

weeks_col = ydata_complete(:,2);

xcal = xdata_complete(ismember(weeks_col,[1,2,3,4,5]),:);
ycal = ydata_complete(ismember(weeks_col,[1,2,3,4,5]),3:end);

xval = xdata_complete(ismember(weeks_col,[6]),:);
yval = ydata_complete(ismember(weeks_col,[6]),3:end);

xtest = xdata_complete(ismember(weeks_col,[7,8]),:);
ytest = ydata_complete(ismember(weeks_col,[7,8]),3:end);


y_labels =  char({'fat', 'lactose', 'protein', 'urea'});

save(mainpath + "/data/d0023_milkrobot_nirsensor2017/data_prepared/trans_milk_nirsensor2017_all_weeks_data_prepared.mat", 'xcal','ycal', 'xval', 'yval','xtest', 'ytest', 'y_labels')



%% individual datasets only calibration per week

for week_ii=1:8
   
    xcal = xdata_complete(ismember(weeks_col,[week_ii]),:);
    ycal = ydata_complete(ismember(weeks_col,[week_ii]),3:end);
    
    save(mainpath + "/data/d0023_milkrobot_nirsensor2017/data_prepared/trans_milk_nirsensor2017_"+"week_" + num2str(week_ii) + "_data_prepared.mat", 'xcal','ycal','y_labels')

end






