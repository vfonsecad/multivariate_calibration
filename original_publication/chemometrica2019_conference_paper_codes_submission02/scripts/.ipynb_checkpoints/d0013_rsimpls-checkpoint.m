% ---------------------------------------------

% ---- model building rsimpls with libra toolbox
% --- by: valeria fonseca diaz

% ---------------------------------------------

clear; clc;

%%
% --- paths and functions
mainpath = "/home/valeria/vfonsecad/kul_phd/conferences/chemometrica-2019/paper/chemometrica2019_conference_paper_codes_submission02";
addpath(mainpath + "/methods/matlab/LIBRA_20160628")

% --- working directories 

% -*-*-*-*-
caseID = "d0013_octane";
matfile = "/octane_data_prepared";

% -*-*-*-*-

wdir_in = mainpath + "/data/"+ caseID + "/data_prepared";
wdir_out = mainpath + "/output";

% --- data

% -*-*-*-*-
iy = 1;
% -*-*-*-*-

load(strcat(wdir_in, matfile, ".mat"));
y_name = y_labels(iy,:) + "_Y" +  string(iy) + "_" + matfile;set(0, 'DefaulttextInterpreter', 'none');
xcal_current = xcal(:,:);
ycal_current = ycal(:,iy);
xval_current = xval(:,:);
yval_current = yval(:,iy);

ncal = size(xcal_current,1);
nval = size(xval_current,1);


disp("data is ready")
%% --- model building cv

% -*-*-*-*-
nlv = 10;
rsimpls_alpha_range = [0.51,0.6,0.7,0.8,0.9];
% -*-*-*-*-


current_rmsecv = zeros(size(rsimpls_alpha_range,2), nlv);
current_r2cv = zeros(size(rsimpls_alpha_range,2), 1);
current_lv_cv = zeros(size(rsimpls_alpha_range,2), 1);
final_rsimpls_alpha_range = zeros(1,size(rsimpls_alpha_range,2));
ii = 1;

for alpha=rsimpls_alpha_range
    
    rsimpls_h = round(alpha*ncal,0) ;
    cv_rsimpls = rsimpls(xcal_current,ycal_current,'rmsecv',1,'plots',0,'kmax',nlv, 'kmaxr', nlv + 1 ,'h', rsimpls_h);
    current_rmsecv(ii,:) = cv_rsimpls.rcs(4,:);
    current_r2cv(ii,:) = cv_rsimpls.rsquared;
    current_lv_cv(ii,:) = cv_rsimpls.k;
    final_rsimpls_alpha_range(1,ii) = cv_rsimpls.h/ncal;
    ii = ii +1;
end
%% --- plot : model building cv 

cmp = colormap(cool(size(current_rmsecv,1)));
for ii=1:size(current_rmsecv,1)
plot(current_rmsecv(ii,:)', 'Color', cmp(ii,:), 'LineWidth', 2)
hold on
end
xticks(1:nlv)
xlabel("lv")
ylabel("error")
ylim([0,1])
title("RSIMPLS - RMSECV", 'FontWeight','normal')
caxis([min(final_rsimpls_alpha_range*100), max(final_rsimpls_alpha_range*100)])
hcb=colorbar;
title(hcb,'% reg. samples')
grid on
set(gca,'FontSize',15)
set(gca,'Position',[0.12 0.14 0.7 0.78]);
saveas(gcf,wdir_out + "/figures/fig_" + caseID + "02_rmsecv_rsimpls.png", 'png')



%% --- selected models

% -*-*-*-*-
selected_models_ids = [1,2,3,4,5];
% -*-*-*-*-


algorithm = cell(length(selected_models_ids),1);
model_id = cell(length(selected_models_ids),1);
rmsep = zeros(length(selected_models_ids),1);
r2p = zeros(length(selected_models_ids),1);
rmsecv = zeros(length(selected_models_ids),1);
r2cv = zeros(length(selected_models_ids),1);
ncp = zeros(length(selected_models_ids),1);
h = zeros(length(selected_models_ids),1);



for ii = 1:length(selected_models_ids)
    
    
    algorithm{ii,1} = "rsimpls";
    model_id{ii,1} = selected_models_ids(ii);
    
    chosen_nlv = current_lv_cv(selected_models_ids(ii));
    chosen_alpha = final_rsimpls_alpha_range(selected_models_ids(ii));

    
    chosen_rsimpls_h = round(chosen_alpha*ncal,0) ;
    mod_rpls = rsimpls(xcal_current,ycal_current,'k',chosen_nlv,'rmsecv',0,'plots',0,'rmsep',1, 'h',chosen_rsimpls_h );
    val_result = predict(xval_current,yval_current, mod_rpls);
    current_rmsep = val_result.rmsep.trim95;
    y_reg = [ones(nval,1) yval_current];
    b = y_reg\val_result.fitted;
    yval_predicted_r2 = y_reg*b;
    current_r2p = 1 - sum((val_result.fitted - yval_predicted_r2).^2)/sum((val_result.fitted - mean(val_result.fitted)).^2);
    
    rmsep(ii) = current_rmsep;
    r2p(ii) = current_r2p;
    rmsecv(ii) = current_rmsecv(selected_models_ids(ii),chosen_nlv);
    r2cv(ii) = current_r2cv(selected_models_ids(ii));
    ncp(ii) = chosen_nlv;
    h(ii) = chosen_alpha;

end
%% selected models summary table

sigma2 = zeros(length(selected_models_ids),1); % only for output display

selected_models_table = table(algorithm,...
model_id,...
ncp,...
h,...
sigma2,...
rmsecv,...
r2cv,...
rmsep,...
r2p);

disp(selected_models_table)


%%  Other results

% --- scores plot


selected_model_id = 4;

chosen_nlv = current_lv_cv(selected_model_id);
chosen_alpha = final_rsimpls_alpha_range(selected_model_id);    
chosen_rsimpls_h = round(chosen_alpha*ncal,0) ;
mod_rpls = rsimpls(xcal_current,ycal_current,'k',chosen_nlv,'rmsecv',0,'plots',0,'rmsep',1, 'h',chosen_rsimpls_h );
   
%%

cmp_id = unique(mod_rpls.flag.all*1);
cmp = colormap(flipud(winter(size(cmp_id,1))));
for ii=1:size(mod_rpls.flag.all,1)
    color_id = find(mod_rpls.flag.all(ii)==cmp_id);
    plot(mod_rpls.T(ii,1),mod_rpls.T(ii,2),'o', 'Color', cmp(color_id,:), 'MarkerSize', 6, 'MarkerFaceColor', cmp(color_id,:))
    hold on
end
xlabel("lv 1")
ylabel("lv 2")
title("RSIMPLS - OCTANE", 'FontWeight', 'normal')
caxis([0,1])
hcb=colorbar;
xlabel(hcb,'sample weights')
grid on
set(gca,'DefaultAxesTitleFontWeight','normal')
set(gca,'FontSize',15)
set(gca,'Position',[0.14 0.14 0.65 0.78]);
saveas(gcf,wdir_out + "/figures/fig_" + caseID + "03_scores_rsimpls_model_0" + num2str(selected_model_id) + ".png")




