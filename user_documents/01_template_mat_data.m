% --------------------------------------------------------------------------

% ----------------- SAVE CHEMOMETRICS DATA  -------------------------------
% by: valeria fonseca diaz
% -------------------------------------------------------------------------

% --- this is a guideline to save the matrices for chemometrics
% --- random matrices are generated
% --- do not change the names of the variables created in the script
% --- only assign existing elements in the workspace as necessary

basepath = "/"; #place here the project directory and change below the corresponding names

mypath = basepath + "/data/my_data/my_data_folder_name/data_prepared";

% --- calibration data (mandatory)

xcal = rand(30,10);
ycal = rand(30,2);

% --- validation data (optional)

xval = rand(22,10);
yval = rand(22,2);


% --- test data (optional)

xtest = rand(36,10);
ytest = rand(36,2);

% --- unlabeled data (optional) 

x_unlabeled = rand(6,10);

% --- y column names (mandatory)

y_labels = char({'y1', 'y2'});

% --- save data (mandatory): exclude from the list below those matrices
% that do not apply (cal and y_names cannot be excluded, exclude from val,
% test or unlabeled

save(mypath + "/mydata.mat", 'y_labels', 'xcal','ycal', 'xval', 'yval', 'xtest','ytest', 'x_unlabeled')

fprintf(" --- SUCCESSFUL --- \n") 
