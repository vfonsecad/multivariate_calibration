% -------------------------------------------------------------
%       create octane data prepared
%by: valeria fonseca diaz
% ------------------------------------------------------------- 
clc ; clear all;

% --- project folder and functions

mainpath = "/home/valeria/vfonsecad/kul_phd/conferences/chemometrica-2019/paper/chemometrica2019_conference_paper_codes_submission02";
addpath(mainpath + "/methods/matlab")


% --- load data

load octane

x = spec.data;
y = octane.data;


% --- select a few samples for validation with kennard stone

KS_val = kennard_stone(x,7);
KS_val.sample_id(find(x(:,200)>0.15)) = 0; % -- keep outliers in the calibration set

% --- get samples into cal and val

xcal = x(KS_val.sample_id==0,:);
ycal = y(KS_val.sample_id==0,:);
xval = x(KS_val.sample_id==1,:);
yval = y(KS_val.sample_id==1,:);

% --- check split

plot(xcal', 'Color', 'blue')
hold on;
plot(xval', 'Color', 'red')
hold on;
hold off;


% --- save data_prepared

y_labels = char({'octane_number'});

save(mainpath + "/data/d0013_octane/data_prepared/octane_data_prepared.mat", 'xcal','ycal', 'xval', 'yval', 'y_labels')


fprintf("finished \n")


