clear all, close all, clc

SET=[];
Time_Sample_Tot=[];
Time_White_Tot=[];
Time_Dark_Tot=[];
Refl_Sample_Tot=[];
Trans_Sample_Tot=[];
Refl_White_Tot=[];
Trans_White_Tot=[];
Refl_Dark_Tot=[];
Trans_Dark_Tot=[];
Refl_Tot=[];
Trans_Tot=[];

Time_Milk_Tot=[];
Cow_ID_Tot=[];
Milk_yield_Tot=[];
Time_PrevMilk_Tot=[];
Milk_Interv_Tot=[];
Fat_Tot=[];
Prot_Tot=[];
Lact_Tot=[];
Urea_Tot=[];
SCC_Tot=[];

pathFolder='C:\Users\bena\Box Sync\ZAPGeel\Research\Datasets\NIRSensor_2017\';
folders=regexp(genpath(pathFolder),['[^;]*'],'match');
for h=1:8
    pathstr = [folders{find(not(cellfun('isempty',strfind(folders,[pathFolder 'SET' num2str(h) '_']))))} '\'];
    % Transform Labview date-number to matlab date and time
    % datestr((Dark1(1,1)/(3600*24))+(2/24)+695422)
    disp(['SET' num2str(h)])
    disp(['-Read files'])
    for i = 1:8
        switch i
            case 1
                str='Dark1';
            case 2
                str='Dark2';
            case 3
                str='White1';
            case 4
                str='White2';
            case 5
                str='Sample1';
            case 6
                str='Sample2';
            case 7
                str='MPR';
            case 8
                str='Robot';
        end
        switch i
            case num2cell(1:6)
                fileID = fopen([pathstr str '.txt']);
                temp=textscan(fileID,repmat('%s ',1,258),'collectoutput',1);
                fclose(fileID);
                %             eval([str '=cell2mat(cellfun(@(x) str2double(strrep(x,'','',''.'')), temp{:},''UniformOutput'',false));']);
                temp2=cell2mat(cellfun(@(x) str2double(strrep(x,',','.')), temp{:},'UniformOutput',false));
                % Transform Labview date-number to matlab date and time (to be used
                % with matlab functions datestr and datenum)
                temp2(:,1)=(temp2(:,1)./(3600.*24))+(2./24)+695422;
                % Remove double lines of specra from the same cow
                temp3=find(diff(temp2(:,2))==0);
                if ~isempty(temp3)
                    warning(['Removed rows [' num2str(temp3') '] from variable ' str]);
                    temp2(temp3,:)=[];
                end
            case {7,8}
                %             eval([str '=xlsread(''' pathstr str '.xlsx'');']);
                temp2=xlsread([pathstr str '.xlsx']);
                if i==7
                    if all(diff(temp2(:,5)))
                        temp2(:,5)=[];
                    else
                        error('Column 5 in ''MPR'' is not an integer series as expected')
                    end
                    if all(all(isnan(temp2(:,2))))
                        temp2(:,2)=[];
                    else
                        error('Column 2 in ''MPR'' has not all NaN-values as expected')
                    end
                    if all((temp2(:,1)-temp2(1,1))==0)
                        temp2(:,1)=[];
                    else
                        error('Column 1 in ''MPR'' has not all the same values as expected')
                    end
                elseif i==8
                    % Transform Excel date-number to matlab date and time (to be used
                    % with matlab functions datestr and datenum)
                    temp2(:,1)=temp2(:,1)+693960;
                    temp2(:,8)=temp2(:,8)+693960;
                    % temp2(:,9) is milking interval in seconds
                    % Length of previous milking in seconds:
                    % hist(((temp2(:,1)-temp2(:,8))-(temp2(:,9)/3600/24))*24*60)
                    if all(all(isnan(temp2(:,[2,7,10]))))
                        temp2(:,[2,7,10])=[];
                    else
                        error('Columns [2, 7, 10] in ''MPR'' do not all have NaN-values as expected')
                    end
                    if all(temp2(:,4))
                        temp2(:,4)=[];
                    else
                        error('Column 4 in ''MPR'' has not all one-values as expected')
                    end
                end
        end
        eval([str '=temp2;']);
    end
    disp(['-Compare cow numbers and length of spectral data'])
    if ~isequal(size(Dark1,1),size(Dark2,1),size(Sample1,1),size(Sample2,1),size(White1,1),size(White2,1))
        error('Spectral files not of same length')
    elseif ~isequal(Dark1(:,2),Dark2(:,2),Sample1(:,2),Sample2(:,2),White1(:,2),White2(:,2))
        error('Spectral files not of the same cows')
    end
    disp(['-Compare cow numbers and length of robot and MPR-data'])
    if size(Robot,1)<size(MPR,1)
        warning(['Number of rows of ''Robot'' (' num2str(size(Robot,1)) ') is smaller than number of rows in ''MPR'' (' num2str(size(MPR,1)) ')'])
    end
    begin=1;
    buildRobot=[];
    buildMPR=[];
    for j=1:size(MPR,1)
        I=min(find(Robot(begin:end,2)==MPR(j,1)))+begin-1;
        if isempty(I)
            warning(['Row [' num2str(j) '] removed from MPR: cow [' num2str(MPR(j,1)) '] not found in ''Robot'' on expected position'])
        else
            if I~=begin
                if ~(mod(Robot(I,4),10)==mod(MPR(j,2),10))
                    error(['Sample number on row [' num2str(I) '] of ''Robot'' is not matching [mod(X,10)] the sample number in row [' num2str(j) '] of ''MPR'''])
                else
                    warning(['Removed rows [' num2str([begin:(I-1)]) '] from ''Robot'''])
                end
            elseif ~(Robot(I,4)==MPR(j,2))
                error(['Sample number on row [' num2str(I) '] of ''Robot'' not the same as sample number in row [' num2str(j) '] of ''MPR'''])
            end
            buildRobot=[buildRobot;I];
            buildMPR=[buildMPR;j];
            begin=I+1;
        end
    end
    if I~=size(Robot,1)
        warning(['Removed rows [' num2str([(I+1):size(Robot,1)]) '] from ''Robot'' (at the end)'])
    end
    Robot=Robot(buildRobot,:);
    MPR=MPR(buildMPR,:);
    figure(),subplot(1,2,1),plot(Robot(:,2)-MPR(:,1))
    title('Diff. cow nr ''Robot'' - ''MPR''')
    subplot(1,2,2),plot(Robot(:,4)-MPR(:,2))
    title('Diff. sample nr ''Robot'' - ''MPR''')
    disp(['-Check milk production levels'])
    LowProduction=find(Robot(:,3)<=6);
    if ~isempty(LowProduction)
        warning(['Removed rows [' num2str(LowProduction') '] from variable ''Robot'' and ''MPR'' because of low milk production']);
        Robot(LowProduction,:)=[];
        MPR(LowProduction,:)=[];
    end
    disp(['-Compare spectral data with robot- and MPR-data'])
    begin=1;
    buildSpectra=[];
    buildRefdata=[];
    for j=1:size(Sample1,1)
        I=find(Robot(begin:end,2)==Sample1(j,2))+begin-1;
        if isempty(I)
            warning(['Row [' num2str(j) '] removed from spectra: cow [' num2str(Sample1(j,2)) '] not found in ''Robot'''])
        else
            [C,J]=min(abs(Robot(I,1)-Sample1(j,1)));
            if C>0.02
                warning(['Row [' num2str(j) '] removed from spectra: no cow [' num2str(Sample1(j,2)) '] in ''Robot'' which matches in a max. time interval of 30 min.'])
            else
                K=I(J); 
                if K~=begin
                    warning(['Removed rows [' num2str([begin:(K-1)]) '] from ''Robot'' and ''MPR'''])
                end
                buildRefdata=[buildRefdata;K];
                buildSpectra=[buildSpectra;j];
                begin=K+1;
            end
        end
    end
    if K~=size(Robot,1)
        warning(['Removed rows [' num2str([(K+1):size(Robot,1)]) '] from ''Robot'' and ''MPR'' (at the end)'])
    end
    disp(['Final checkup'])
    if ~isequal(length(buildSpectra),length(buildRefdata))
        error('Selection vectors not of same length')
    elseif ~isequal(Sample1(buildSpectra,2),Robot(buildRefdata,2),MPR(buildRefdata,1))
        error('Files do not contain the same cows')
    end

    if h==2
        % The rack with samples 111 - 120 (cow 7679 - cow 1867) was likely
        % switched with the rack with samples 121 - 130 (cow 3803 - cow
        % 3808). Switching the reference values results in much better
        % prediction of the model
        temp = MPR(111:120,3:7);
        MPR(111:120,3:7)=MPR(121:130,3:7);
        MPR(121:130,3:7)=temp;
    elseif h==4
        % Cuvette was not properly filled for sample 23 of SET 4 (=Cow
        % 7488) resulting in low Reflectance and high Transmittance
        % (saturation) of sample spectrum. This sample is in position 9
        % after previous filtering.
        buildRefdata(9)=[];
        buildSpectra(9)=[];
    elseif h==7
        % White spectra of both Reflectance and Transmittance were
        % significantly different for sample 116 of SET 7 (= Cow 910) for unknown
        % reasons. This sample is in position 94
        % after previous filtering.
        buildRefdata(94)=[];
        buildSpectra(94)=[];
    end
    
    disp(['---------------------'])
    Time_Sample=Sample1(buildSpectra,1);
    Time_White=White1(buildSpectra,1);
    Time_Dark=Dark1(buildSpectra,1);
    Refl_Sample=Sample1(buildSpectra,3:end);
    Trans_Sample=Sample2(buildSpectra,3:end);
    Refl_White=White1(buildSpectra,3:end);
    Trans_White=White2(buildSpectra,3:end);
    Refl_Dark=Dark1(buildSpectra,3:end);
    Trans_Dark=Dark2(buildSpectra,3:end);
    Refl=((Refl_Sample-repmat(mean(Refl_Dark),size(Refl_Dark,1),1))./(repmat(mean(Refl_White),size(Refl_White,1),1)-repmat(mean(Refl_Dark),size(Refl_Dark,1),1)))+0.015;
    Trans=(Trans_Sample-repmat(mean(Trans_Dark),size(Trans_Dark,1),1))./(repmat(mean(Trans_White),size(Trans_White,1),1)-repmat(mean(Trans_Dark),size(Trans_Dark,1),1));
    
    Time_Milk=Robot(buildRefdata,1);
    Cow_ID=Robot(buildRefdata,2);
    Milk_yield=Robot(buildRefdata,3);
    Time_PrevMilk=Robot(buildRefdata,5);
    Milk_Interv=Robot(buildRefdata,6); % in seconds!!
    Fat=MPR(buildRefdata,3);
    Prot=MPR(buildRefdata,4);
    Lact=MPR(buildRefdata,5);
    Urea=MPR(buildRefdata,6);
    SCC=MPR(buildRefdata,7);
    
    save([pathstr 'Data_SET' num2str(h) '.mat'],'pathFolder','Time_Sample','Time_White','Time_Dark','Refl_Sample','Trans_Sample','Refl_White','Trans_White','Refl_Dark','Trans_Dark','Refl','Trans','Time_Milk','Cow_ID','Milk_yield','Time_PrevMilk','Milk_Interv','Fat','Prot','Lact','Urea','SCC');
    
    SET=[SET;repmat(h,length(buildRefdata),1)];
    Time_Sample_Tot=[Time_Sample_Tot;Time_Sample];
    Time_White_Tot=[Time_White_Tot;Time_White];
    Time_Dark_Tot=[Time_Dark_Tot;Time_Dark];
    Refl_Sample_Tot=[Refl_Sample_Tot;Refl_Sample];
    Trans_Sample_Tot=[Trans_Sample_Tot;Trans_Sample];
    Refl_White_Tot=[Refl_White_Tot;Refl_White];
    Trans_White_Tot=[Trans_White_Tot;Trans_White];
    Refl_Dark_Tot=[Refl_Dark_Tot;Refl_Dark];
    Trans_Dark_Tot=[Trans_Dark_Tot;Trans_Dark];
    Refl_Tot=[Refl_Tot;Refl];
    Trans_Tot=[Trans_Tot;Trans];
    
    Time_Milk_Tot=[Time_Milk_Tot;Time_Milk];
    Cow_ID_Tot=[Cow_ID_Tot;Cow_ID];
    Milk_yield_Tot=[Milk_yield_Tot;Milk_yield];
    Time_PrevMilk_Tot=[Time_PrevMilk_Tot;Time_PrevMilk];
    Milk_Interv_Tot=[Milk_Interv_Tot;Milk_Interv];
    Fat_Tot=[Fat_Tot;Fat];
    Prot_Tot=[Prot_Tot;Prot];
    Lact_Tot=[Lact_Tot;Lact];
    Urea_Tot=[Urea_Tot;Urea];
    SCC_Tot=[SCC_Tot;SCC];
    clearvars -except pathFolder folders h SET Time_Sample_Tot Time_White_Tot Time_Dark_Tot Refl_Sample_Tot Trans_Sample_Tot Refl_White_Tot Trans_White_Tot Refl_Dark_Tot Trans_Dark_Tot Refl_Tot Trans_Tot Time_Milk_Tot Cow_ID_Tot Milk_yield_Tot Time_PrevMilk_Tot Milk_Interv_Tot Fat_Tot Prot_Tot Lact_Tot Urea_Tot SCC_Tot
end
save([pathFolder 'Data_Tot.mat'],'pathFolder','SET','Time_Sample_Tot','Time_White_Tot','Time_Dark_Tot','Refl_Sample_Tot','Trans_Sample_Tot','Refl_White_Tot','Trans_White_Tot','Refl_Dark_Tot','Trans_Dark_Tot','Refl_Tot','Trans_Tot','Time_Milk_Tot','Cow_ID_Tot','Milk_yield_Tot','Time_PrevMilk_Tot','Milk_Interv_Tot','Fat_Tot','Prot_Tot','Lact_Tot','Urea_Tot','SCC_Tot');
% figure(),plot(isnan(Fat_Tot)+isnan(Prot_Tot)+isnan(Lact_Tot)+isnan(Urea_Tot)+isnan(SCC_Tot))

%%
figure()
for i=1:6
    switch i
        case 1
            plotdata=Refl_Dark_Tot;
            titleline='Dark - Reflectance';
        case 2
            plotdata=Refl_White_Tot;
            titleline='White - Reflectance';
        case 3
            plotdata=Refl_Sample_Tot;
            titleline='Sample - Reflectance';
        case 4            
            plotdata=Trans_Dark_Tot;
            titleline='Dark - Transmittance';
        case 5
            plotdata=Trans_White_Tot;
            titleline='White - Transmittance';
        case 6
            plotdata=Trans_Sample_Tot;
            titleline='Sample - Transmittance';
    end
    subplot(2,3,i)
    imagesc((plotdata./repmat(mean(plotdata),size(plotdata,1),1)))
    [C,ia,ic]=unique(SET);
%     set(gca,'YTick',ia,'YTickLabel',C);
%     set(gca,'XTick',[24:33:256],'XTickLabel',[1000:100:1700]);
    xlabel('Wavelengths')
    ylabel('SET')
    title(titleline)
end

figure()
for i=1:6
    switch i
        case 1
            plotdata=Refl_Dark_Tot;
            titleline='Dark - Reflectance';
        case 2
            plotdata=Refl_White_Tot;
            titleline='White - Reflectance';
        case 3
            plotdata=Refl_Sample_Tot;
            titleline='Sample - Reflectance';
        case 4            
            plotdata=Trans_Dark_Tot;
            titleline='Dark - Transmittance';
        case 5
            plotdata=Trans_White_Tot;
            titleline='White - Transmittance';
        case 6
            plotdata=Trans_Sample_Tot;
            titleline='Sample - Transmittance';
    end
    subplot(2,3,i)
    plot(plotdata')
    set(gca,'XTick',[24:33:256],'XTickLabel',[1000:100:1700]);
    xlabel('Wavelengths')
    ylabel('Value')
    title(titleline)
end

%%
Dark1Av=mean(Dark1,1);
Dark2Av=mean(Dark2,1);
White1Av=mean(White1,1);
White2Av=mean(White2,1);

% figure(),plot(White1(:,120))
% figure(),surf(White1(:,3:end)./repmat(White1Av(:,3:end),size(White1,1),1));

temp=Dark1(:,3:end)./repmat(Dark1Av(:,3:end),size(Dark1,1),1);
figure(),plot([mean(temp(:,20:235),2),savgol(mean(temp(:,20:235),2)',15)']);
Dark1Av2=mean(temp(:,20:235),2)*Dark1Av;
Dark1Av3=savgol(mean(temp(:,20:235),2)',15)'*Dark1Av;

temp=Dark2(:,3:end)./repmat(Dark2Av(:,3:end),size(Dark2,1),1);
figure(),plot([mean(temp(:,20:235),2),savgol(mean(temp(:,20:235),2)',15)']);
Dark2Av2=mean(temp(:,20:235),2)*Dark2Av;
Dark2Av3=savgol(mean(temp(:,20:235),2)',15)'*Dark2Av;

temp=White1(:,3:end)./repmat(White1Av(:,3:end),size(White1,1),1);
figure(),plot([mean(temp(:,20:235),2),savgol(mean(temp(:,20:235),2)',15)']);
White1Av2=mean(temp(:,20:235),2)*White1Av;
White1Av3=savgol(mean(temp(:,20:235),2)',15)'*White1Av;

temp=White2(:,3:end)./repmat(White2Av(:,3:end),size(White2,1),1);
figure(),plot([mean(temp(:,20:235),2),savgol(mean(temp(:,20:235),2)',15)']);
White2Av2=mean(temp(:,20:235),2)*White2Av;
White2Av3=savgol(mean(temp(:,20:235),2)',15)'*White2Av;

Refl=(Sample1(:,3:end)-Dark1(:,3:end))./(White1(:,3:end)-Dark1(:,3:end))+0.015;
Trans=(Sample2(:,3:end)-Dark2(:,3:end))./(White2(:,3:end)-Dark2(:,3:end));

ReflAv=(Sample1(:,3:end)-repmat(Dark1Av(:,3:end),size(Dark1,1),1))./(repmat(White1Av(:,3:end),size(White1,1),1)-repmat(Dark1Av(:,3:end),size(Dark1,1),1))+0.015;
TransAv=(Sample2(:,3:end)-repmat(Dark2Av(:,3:end),size(Dark2,1),1))./(repmat(White2Av(:,3:end),size(White2,1),1)-repmat(Dark2Av(:,3:end),size(Dark2,1),1));

ReflAv2=(Sample1(:,3:end)-Dark1Av2(:,3:end))./(White1Av2(:,3:end)-Dark1Av2(:,3:end))+0.015;
TransAv2=(Sample2(:,3:end)-Dark2Av2(:,3:end))./(White2Av2(:,3:end)-Dark2Av2(:,3:end));

ReflAv3=(Sample1(:,3:end)-Dark1Av3(:,3:end))./(White1Av3(:,3:end)-Dark1Av3(:,3:end))+0.015;
TransAv3=(Sample2(:,3:end)-Dark2Av3(:,3:end))./(White2Av3(:,3:end)-Dark2Av3(:,3:end));

figure(),plot([925:3:1690],Refl')
title('Reflectance')
xlabel('Wavelength (nm)')
ylabel('Reflectance (*100%)')
axis([950,1700,0,0.35])

figure(),plot([925:3:1690],Trans'./100)
title('Transmittance')
xlabel('Wavelength (nm)')
ylabel('Transmittance (*100%)')
axis([950,1700,0,0.08])

%% Test Model_Dries on data
% Transmittance
wavelenghtselection = [3:145 192:258]-2;
path_model=[pathFolder 'Kalibratiemodel\Transmittance\'];

% Fat Model: No transformation to absorbance
g=4;
load([path_model 'Vet_Model.mat'])
Xttrvar = dataset(Trans_Tot(:,wavelenghtselection));
Xttrvar.axisscale{2} = 1:length(wavelenghtselection);
Xttrvar.include{2}=Model_par_use{g};
validvarpls = pls(Xttrvar,zeros(size(Trans_Tot,1),1),Model_par_modelvarpls{g},Model_par_optionspls{g});
Fat_Trans = cell2mat(validvarpls.pred)';

% Protein Model: No transformation to absorbance
g = 4;
load([path_model 'Eiwit_Model.mat'])
Xttrvar = dataset(Trans_Tot(:,wavelenghtselection ));
Xttrvar.axisscale{2} = 1:length(wavelenghtselection);
Xttrvar.include{2}=Model_par_use{g};
validvarpls = pls(Xttrvar,zeros(size(Trans_Tot,1),1),Model_par_modelvarpls{g},Model_par_optionspls{g});
Prot_Trans = cell2mat(validvarpls.pred)';

% Lactose Model: Transformation to absorbance
g = 4;
load([path_model 'Lactose_Model.mat'])
Xttrvar = dataset(log10(1./Trans_Tot(:,wavelenghtselection )));
Xttrvar.axisscale{2} = 1:length(wavelenghtselection);
Xttrvar.include{2}=Model_par_use{g};
validvarpls = pls(Xttrvar,zeros(size(Trans_Tot,1),1),Model_par_modelvarpls{g},Model_par_optionspls{g});
Lact_Trans  = cell2mat(validvarpls.pred)';

figure()
subplot(2,3,1),plot(Fat_Tot,Fat_Trans','go',[0,10],[0,10],'k-')
legend({['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((Fat_Trans'-Fat_Tot).^2)))]},'location','NorthWest');
title('Transmittance')
xlabel('Fat (%w/w) by reference')
ylabel('Fat (%w/w) by sensor')
axis([0,8,0,8])

subplot(2,3,2),plot(Prot_Tot,Prot_Trans','go',[0,10],[0,10],'k-')
legend({['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((Prot_Trans'-Prot_Tot).^2)))]},'location','NorthWest');
title('Transmittance')
xlabel('Protein (%w/w) by reference')
ylabel('Protein (%w/w) by sensor')
axis([2.5,4.5,1,3.5])

subplot(2,3,3),plot(Lact_Tot,Lact_Trans','go',[0,10],[0,10],'k-')
legend({['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((Lact_Trans'-Lact_Tot).^2)))]},'location','NorthWest');
title('Transmittance')
xlabel('Lactose (%w/w) by reference')
ylabel('Lactose (%w/w) by sensor')
axis([3.5,5.5,-3,3])

% Reflectance
wavelenghtselection = [3:145 192:258]-2;
path_model=[pathFolder 'Kalibratiemodel\Reflectance\'];

% Fat Model: Transformation to absorbance
g=5;
load([path_model 'Vet_Model.mat'])
Xttrvar = dataset(log10(1./Refl_Tot(:,wavelenghtselection )));
Xttrvar.axisscale{2} = 1:length(wavelenghtselection);
Xttrvar.include{2}=Model_par_use{g};
validvarpls = pls(Xttrvar,zeros(size(Refl_Tot,1),1),Model_par_modelvarpls{g},Model_par_optionspls{g});
Fat_Refl = cell2mat(validvarpls.pred)';

% Protein Model: Transformation to absorbance
g = 5;
load([path_model 'Eiwit_Model.mat'])
Xttrvar = dataset(log10(1./Refl_Tot(:,wavelenghtselection )));
Xttrvar.axisscale{2} = 1:length(wavelenghtselection);
Xttrvar.include{2}=Model_par_use{g};
validvarpls = pls(Xttrvar,zeros(size(Refl_Tot,1),1),Model_par_modelvarpls{g},Model_par_optionspls{g});
Prot_Refl = cell2mat(validvarpls.pred)';

% Lactose Model: Transformation to absorbance
g = 4;
load([path_model 'Lactose_Model.mat'])
Xttrvar = dataset(log10(1./Refl_Tot(:,wavelenghtselection )));
Xttrvar.axisscale{2} = 1:length(wavelenghtselection);
Xttrvar.include{2}=Model_par_use{g};
validvarpls = pls(Xttrvar,zeros(size(Refl_Tot,1),1),Model_par_modelvarpls{g},Model_par_optionspls{g});
Lact_Refl  = cell2mat(validvarpls.pred)';

subplot(2,3,4),plot(Fat_Tot,Fat_Refl','go',[-2,10],[-2,10],'k-')
legend({['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((Fat_Refl'-Fat_Tot).^2)))]},'location','NorthWest');
title('Transmittance')
xlabel('Fat (%w/w) by reference')
ylabel('Fat (%w/w) by sensor')
axis([1,8,-3,6])

subplot(2,3,5),plot(Prot_Tot,Prot_Refl','go',[0,10],[0,10],'k-')
legend({['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((Prot_Refl'-Prot_Tot).^2)))]},'location','NorthWest');
title('Reflectance')
xlabel('Protein (%w/w) by reference')
ylabel('Protein (%w/w) by sensor')
axis([2.5,4.5,3,7])

subplot(2,3,6),plot(Lact_Tot,Lact_Refl','go',[0,10],[0,10],'k-')
legend({['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((Lact_Refl'-Lact_Tot).^2)))]},'location','NorthWest');
title('Reflectance')
xlabel('Lactose (%w/w) by reference')
ylabel('Lactose (%w/w) by sensor')
axis([3.5,5.5,3.5,8.5])

%% Same plots as previous, but with different colors for different SET's
figure()
for i=1:6
    subplot(2,3,i)
    plot([0,10],[0,10],'k-')
    switch i
        case 1
            xlabel('Fat (%w/w) by reference')
            ylabel('Fat (%w/w) by sensor')
            title('Transmittance')
            lab_values=Fat_Tot(~isnan(Fat_Tot));
            sens_values=Fat_Trans(~isnan(Fat_Tot))';
            SET_values=SET(~isnan(Fat_Tot));
            axis([1,8,0,8])
        case 2
            xlabel('Protein (%w/w) by reference')
            ylabel('Protein (%w/w) by sensor')
            title('Transmittance')
            lab_values=Prot_Tot(~isnan(Prot_Tot));
            sens_values=Prot_Trans(~isnan(Prot_Tot))';
            SET_values=SET(~isnan(Prot_Tot));
            axis([2.5,4.5,1,3.5])
        case 3
            xlabel('Lactose (%w/w) by reference')
            ylabel('Lactose (%w/w) by sensor')
            title('Transmittance')
            lab_values=Lact_Tot(~isnan(Lact_Tot));
            sens_values=Lact_Trans(~isnan(Lact_Tot))';
            axis([3.5,5.5,-3,3])
            SET_values=SET(~isnan(Lact_Tot));
        case 4
            xlabel('Fat (%w/w) by reference')
            ylabel('Fat (%w/w) by sensor')
            title('Reflectance')
            lab_values=Fat_Tot(~isnan(Fat_Tot));
            sens_values=Fat_Refl(~isnan(Fat_Tot))';
            SET_values=SET(~isnan(Fat_Tot));
            axis([1,8,-3,6])
        case 5
            xlabel('Protein (%w/w) by reference')
            ylabel('Protein (%w/w) by sensor')
            title('Reflectance')
            lab_values=Prot_Tot(~isnan(Prot_Tot));
            sens_values=Prot_Refl(~isnan(Prot_Tot))';
            SET_values=SET(~isnan(Prot_Tot));
            axis([2.5,4.5,3,7])
        case 6
            xlabel('Lactose (%w/w) by reference')
            ylabel('Lactose (%w/w) by sensor')
            title('Reflectance')
            lab_values=Lact_Tot(~isnan(Lact_Tot));
            sens_values=Lact_Refl(~isnan(Lact_Tot))';
            SET_values=SET(~isnan(Lact_Tot));
            axis([3.5,5.5,3.5,8.5])
    end
    hold all
    dot_color={'g*','b*','r*','c*','k*','m*','y*','go'};
    for j=1:8
        plot(lab_values(SET_values==j),sens_values(SET_values==j),dot_color{j})
    end
    hold off
end

%% Same plots as previous, but with bias correction
figure()
for i=1:6
    subplot(2,3,i)
    plot([0,10],[0,10],'k-')
    switch i
        case 1
            xlabel('Fat (%w/w) by reference')
            ylabel('Fat (%w/w) by sensor')
            title('Transmittance')
            lab_values=Fat_Tot(~isnan(Fat_Tot));
            sens_values=Fat_Trans(~isnan(Fat_Tot))';
            SET_values=SET(~isnan(Fat_Tot));
            axis([1,8,0,7])
        case 2
            xlabel('Protein (%w/w) by reference')
            ylabel('Protein (%w/w) by sensor')
            title('Transmittance')
            lab_values=Prot_Tot(~isnan(Prot_Tot));
            sens_values=Prot_Trans(~isnan(Prot_Tot))';
            SET_values=SET(~isnan(Prot_Tot));
            axis([2.5,4.5,2,4.5])
        case 3
            xlabel('Lactose (%w/w) by reference')
            ylabel('Lactose (%w/w) by sensor')
            title('Transmittance')
            lab_values=Lact_Tot(~isnan(Lact_Tot));
            sens_values=Lact_Trans(~isnan(Lact_Tot))';
            axis([3.5,5.5,1.5,7.5])
            SET_values=SET(~isnan(Lact_Tot));
        case 4
            xlabel('Fat (%w/w) by reference')
            ylabel('Fat (%w/w) by sensor')
            title('Reflectance')
            lab_values=Fat_Tot(~isnan(Fat_Tot));
            sens_values=Fat_Refl(~isnan(Fat_Tot))';
            SET_values=SET(~isnan(Fat_Tot));
            axis([1,8,0,8])
        case 5
            xlabel('Protein (%w/w) by reference')
            ylabel('Protein (%w/w) by sensor')
            title('Reflectance')
            lab_values=Prot_Tot(~isnan(Prot_Tot));
            sens_values=Prot_Refl(~isnan(Prot_Tot))';
            SET_values=SET(~isnan(Prot_Tot));
            axis([2.5,4.5,1,5])
        case 6
            xlabel('Lactose (%w/w) by reference')
            ylabel('Lactose (%w/w) by sensor')
            title('Reflectance')
            lab_values=Lact_Tot(~isnan(Lact_Tot));
            sens_values=Lact_Refl(~isnan(Lact_Tot))';
            SET_values=SET(~isnan(Lact_Tot));
            axis([3.5,5.5,2.5,7])
    end
    hold all
    p = polyfit(lab_values,sens_values,1);
    sens_values2=sens_values-(polyval(p,lab_values)-lab_values); 
    dot_color={'g*','b*','r*','c*','k*','m*','y*','go'};
    legend_values{1}=['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((sens_values2-lab_values).^2)))];
    for j=1:8
        legend_values{j+1}=['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((sens_values2(SET_values==j)-lab_values(SET_values==j)).^2)))];
        plot(lab_values(SET_values==j),sens_values2(SET_values==j),dot_color{j})
    end
    hold off
    legend(legend_values,'location','NorthWest');
end

%% Recalculate PLS-model with same pre-processing and variable selection...
% Transmittance
ncomp=20;
wavelenghtselection = [3:145 192:258]-2;
path_model=[pathFolder 'Kalibratiemodel\Transmittance\'];
[~,ia,~]=unique(SET);
I=repmat((ia(2:end))',2,1);
J=repmat([-10 10 10 -10],1,length(ia));

% Fat Model: No transformation to absorbance
g=4;
load([path_model 'Vet_Model.mat'])
Y=Fat_Tot(1:150);
Yt=Fat_Tot(151:end);
Xttrvar = dataset(Trans_Tot(:,wavelenghtselection));
Xttrvar.axisscale{2} = 1:length(wavelenghtselection);
Xttrvar.include{2}=Model_par_use{g};
[press,cumpress,rmsecv,rmsec,cvpred,misclassed] = crossval(Xttrvar(1:150,:),Y,'sim',Cow_ID_Tot(1:150),ncomp,Model_par_optionspls{g});
% Calculate absolute residuals
for j = 1:ncomp
    absres(:,j) = abs(cvpred(:,j)-Y);
end
% select number of latent variables based on F test
[cumpressmin,h_star] = min(cumpress);
Hlv = 1;
h = 1;
F_amm = ftest(0.25,length(Y),length(Y));
while Hlv == 1;
    F = cumpress(h)/cumpress(h_star);
    if F < F_amm
        Hlv = 0;
    else
        h = h+1;
    end
end
nsel = h;
clear h Hlv F h_star F_amm;
modelvarpls = pls(Xttrvar(1:150,:),Y,nsel,Model_par_optionspls{g});
validvarpls = pls(Xttrvar(151:end,:),Yt,modelvarpls,Model_par_optionspls{g});
figure(99),subplot(2,3,1),plot(Y,cvpred(:,nsel)','b*',Yt,cell2mat(validvarpls.pred)','go',[0,10],[0,10],'k-')
legend({['RMSE_C_V = ' sprintf('%0.3f',rmsecv(nsel))],['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((cell2mat(validvarpls.pred)-Yt).^2)))]},'location','NorthWest','FontSize',12);
title(['Transmittance - LV = ' num2str(nsel)],'FontSize',14)
xlabel('Fat (%w/w) by reference','FontSize',14)
ylabel('Fat (%w/w) by sensor','FontSize',14)
axis([1,8,1,8])
figure(100),subplot(2,3,1),plot([1:length(cvpred(:,nsel))],cvpred(:,nsel)-Y,[1:length(Yt)]+length(cvpred(:,nsel)),(cell2mat(validvarpls.pred)-Yt),[0,1300],[0,0],'k-')
line(I(:)',J(1:length(I(:))),'color','k','LineWidth',2);
ylabel('Error (%w/w)    [Y_{Pred} - Y_{Ref}]','FontSize',14)
title('Transmittance - Fat','FontSize',14)
xlabel('Samples','FontSize',14)
ylim([-1,1])
xlim([0,1300])

% Protein Model: No transformation to absorbance
g = 4;
load([path_model 'Eiwit_Model.mat'])
Y=Prot_Tot(1:150);
Yt=Prot_Tot(151:end);
Xttrvar = dataset(Trans_Tot(:,wavelenghtselection ));
Xttrvar.axisscale{2} = 1:length(wavelenghtselection);
Xttrvar.include{2}=Model_par_use{g};
[press,cumpress,rmsecv,rmsec,cvpred,misclassed] = crossval(Xttrvar(1:150,:),Y,'sim',Cow_ID_Tot(1:150),ncomp,Model_par_optionspls{g});
% Calculate absolute residuals
for j = 1:ncomp
    absres(:,j) = abs(cvpred(:,j)-Y);
end
% select number of latent variables based on F test
[cumpressmin,h_star] = min(cumpress);
Hlv = 1;
h = 1;
F_amm = ftest(0.25,length(Y),length(Y));
while Hlv == 1;
    F = cumpress(h)/cumpress(h_star);
    if F < F_amm
        Hlv = 0;
    else
        h = h+1;
    end
end
nsel = h;
clear h Hlv F h_star F_amm;
modelvarpls = pls(Xttrvar(1:150,:),Y,nsel,Model_par_optionspls{g});
validvarpls = pls(Xttrvar(151:end,:),Yt,modelvarpls,Model_par_optionspls{g});
Eiwit_Trans2 = cvpred(:,nsel);
Eiwit_Trans2t = cell2mat(validvarpls.pred)';
figure(99),subplot(2,3,2),plot(Y,cvpred(:,nsel)','b*',Yt,cell2mat(validvarpls.pred)','go',[0,10],[0,10],'k-')
legend({['RMSE_C_V = ' sprintf('%0.3f',rmsecv(nsel))],['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((cell2mat(validvarpls.pred)-Yt).^2)))]},'location','NorthWest','FontSize',12);
title(['Transmittance - LV = ' num2str(nsel)],'FontSize',14)
xlabel('Protein (%w/w) by reference','FontSize',14)
ylabel('Protein (%w/w) by sensor','FontSize',14)
axis([2,5,2,5])
figure(100),subplot(2,3,2),plot([1:length(cvpred(:,nsel))],cvpred(:,nsel)-Y,[1:length(Yt)]+length(cvpred(:,nsel)),(cell2mat(validvarpls.pred)-Yt),[0,1300],[0,0],'k-')
line(I(:)',J(1:length(I(:))),'color','k','LineWidth',2);
ylabel('Error (%w/w)    [Y_{Pred} - Y_{Ref}]','FontSize',14)
title('Transmittance - Protein','FontSize',14)
xlabel('Samples','FontSize',14)
ylim([-1,1])
xlim([0,1300])

% Lactose Model: Transformation to absorbance
g = 4;
load([path_model 'Lactose_Model.mat'])
Y=Lact_Tot(1:150);
Yt=Lact_Tot(151:end);
Xttrvar = dataset(log10(1./Trans_Tot(:,wavelenghtselection )));
Xttrvar.axisscale{2} = 1:length(wavelenghtselection);
Xttrvar.include{2}=Model_par_use{g};
[press,cumpress,rmsecv,rmsec,cvpred,misclassed] = crossval(Xttrvar(1:150,:),Y,'sim',Cow_ID_Tot(1:150),ncomp,Model_par_optionspls{g});
% Calculate absolute residuals
for j = 1:ncomp
    absres(:,j) = abs(cvpred(:,j)-Y);
end
% select number of latent variables based on F test
[cumpressmin,h_star] = min(cumpress);
Hlv = 1;
h = 1;
F_amm = ftest(0.25,length(Y),length(Y));
while Hlv == 1;
    F = cumpress(h)/cumpress(h_star);
    if F < F_amm
        Hlv = 0;
    else
        h = h+1;
    end
end
nsel = h;
clear h Hlv F h_star F_amm;
modelvarpls = pls(Xttrvar(1:150,:),Y,nsel,Model_par_optionspls{g});
validvarpls = pls(Xttrvar(151:end,:),Yt,modelvarpls,Model_par_optionspls{g});
figure(99),subplot(2,3,3),plot(Y,cvpred(:,nsel)','b*',Yt,cell2mat(validvarpls.pred)','go',[0,10],[0,10],'k-')
legend({['RMSE_C_V = ' sprintf('%0.3f',rmsecv(nsel))],['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((cell2mat(validvarpls.pred)-Yt).^2)))]},'location','NorthWest','FontSize',12);
title(['Transmittance - LV = ' num2str(nsel)],'FontSize',14)
xlabel('Lactose (%w/w) by reference','FontSize',14)
ylabel('Lactose (%w/w) by sensor','FontSize',14)
axis([3.5,5.5,3.5,5.5])
figure(100),subplot(2,3,3),plot([1:length(cvpred(:,nsel))],cvpred(:,nsel)-Y,[1:length(Yt)]+length(cvpred(:,nsel)),(cell2mat(validvarpls.pred)-Yt),[0,1300],[0,0],'k-')
line(I(:)',J(1:length(I(:))),'color','k','LineWidth',2);
ylabel('Error (%w/w)    [Y_{Pred} - Y_{Ref}]','FontSize',14)
title('Transmittance - Lactose','FontSize',14)
xlabel('Samples','FontSize',14)
ylim([-1,1])
xlim([0,1300])
    
% Reflectance
wavelenghtselection = [3:145 192:258]-2;
path_model=[pathFolder 'Kalibratiemodel\Reflectance\'];

% Fat Model: Transformation to absorbance
g=5;
load([path_model 'Vet_Model.mat'])
Y=Fat_Tot(1:150);
Yt=Fat_Tot(151:end);
Xttrvar = dataset(log10(1./Refl_Tot(:,wavelenghtselection )));
Xttrvar.axisscale{2} = 1:length(wavelenghtselection);
Xttrvar.include{2}=Model_par_use{g};
[press,cumpress,rmsecv,rmsec,cvpred,misclassed] = crossval(Xttrvar(1:150,:),Y,'sim',Cow_ID_Tot(1:150),ncomp,Model_par_optionspls{g});
% Calculate absolute residuals
for j = 1:ncomp
    absres(:,j) = abs(cvpred(:,j)-Y);
end
% select number of latent variables based on F test
[cumpressmin,h_star] = min(cumpress);
Hlv = 1;
h = 1;
F_amm = ftest(0.25,length(Y),length(Y));
while Hlv == 1;
    F = cumpress(h)/cumpress(h_star);
    if F < F_amm
        Hlv = 0;
    else
        h = h+1;
    end
end
nsel = h;
clear h Hlv F h_star F_amm;
modelvarpls = pls(Xttrvar(1:150,:),Y,nsel,Model_par_optionspls{g});
validvarpls = pls(Xttrvar(151:end,:),Yt,modelvarpls,Model_par_optionspls{g});
Vet_Refl2 = cvpred(:,nsel);
Vet_Refl2t = cell2mat(validvarpls.pred)';
figure(99),subplot(2,3,4),plot(Y,cvpred(:,nsel)','b*',Yt,cell2mat(validvarpls.pred)','go',[0,10],[0,10],'k-')
legend({['RMSE_C_V = ' sprintf('%0.3f',rmsecv(nsel))],['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((cell2mat(validvarpls.pred)-Yt).^2)))]},'location','NorthWest','FontSize',12);
title(['Reflectance - LV = ' num2str(nsel)],'FontSize',14)
xlabel('Fat (%w/w) by reference','FontSize',14)
ylabel('Fat (%w/w) by sensor','FontSize',14)
axis([1,8,1,8])
figure(100),subplot(2,3,4),plot([1:length(cvpred(:,nsel))],cvpred(:,nsel)-Y,[1:length(Yt)]+length(cvpred(:,nsel)),(cell2mat(validvarpls.pred)-Yt),[0,1300],[0,0],'k-')
line(I(:)',J(1:length(I(:))),'color','k','LineWidth',2);
ylabel('Error (%w/w)    [Y_{Pred} - Y_{Ref}]','FontSize',14)
title('Reflectance - Fat','FontSize',14)
xlabel('Samples','FontSize',14)
ylim([-1,1])
xlim([0,1300])

% Protein Model: Transformation to absorbance
g = 5;
load([path_model 'Eiwit_Model.mat'])
Y=Prot_Tot(1:150);
Yt=Prot_Tot(151:end);
Xttrvar = dataset(log10(1./Refl_Tot(:,wavelenghtselection )));
Xttrvar.axisscale{2} = 1:length(wavelenghtselection);
Xttrvar.include{2}=Model_par_use{g};
[press,cumpress,rmsecv,rmsec,cvpred,misclassed] = crossval(Xttrvar(1:150,:),Y,'sim',Cow_ID_Tot(1:150),ncomp,Model_par_optionspls{g});
% Calculate absolute residuals
for j = 1:ncomp
    absres(:,j) = abs(cvpred(:,j)-Y);
end
% select number of latent variables based on F test
[cumpressmin,h_star] = min(cumpress);
Hlv = 1;
h = 1;
F_amm = ftest(0.25,length(Y),length(Y));
while Hlv == 1;
    F = cumpress(h)/cumpress(h_star);
    if F < F_amm
        Hlv = 0;
    else
        h = h+1;
    end
end
nsel = h;
clear h Hlv F h_star F_amm;
modelvarpls = pls(Xttrvar(1:150,:),Y,nsel,Model_par_optionspls{g});
validvarpls = pls(Xttrvar(151:end,:),Yt,modelvarpls,Model_par_optionspls{g});
Eiwit_Refl2 = cvpred(:,nsel);
Eiwit_Refl2t = cell2mat(validvarpls.pred)';
figure(99),subplot(2,3,5),plot(Y,cvpred(:,nsel)','b*',Yt,cell2mat(validvarpls.pred)','go',[0,10],[0,10],'k-')
legend({['RMSE_C_V = ' sprintf('%0.3f',rmsecv(nsel))],['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((cell2mat(validvarpls.pred)-Yt).^2)))]},'location','NorthWest','FontSize',12);
title(['Reflectance - LV = ' num2str(nsel)],'FontSize',14)
xlabel('Protein (%w/w) by reference','FontSize',14)
ylabel('Protein (%w/w) by sensor','FontSize',14)
axis([2,5,2,5])
figure(100),subplot(2,3,5),plot([1:length(cvpred(:,nsel))],cvpred(:,nsel)-Y,[1:length(Yt)]+length(cvpred(:,nsel)),(cell2mat(validvarpls.pred)-Yt),[0,1300],[0,0],'k-')
line(I(:)',J(1:length(I(:))),'color','k','LineWidth',2);
ylabel('Error (%w/w)    [Y_{Pred} - Y_{Ref}]','FontSize',14)
title('Reflectance - Protein','FontSize',14)
xlabel('Samples','FontSize',14)
ylim([-1,1])
xlim([0,1300])

% Lactose Model: Transformation to absorbance
g = 4;
load([path_model 'Lactose_Model.mat'])
Y=Lact_Tot(1:150);
Yt=Lact_Tot(151:end);
Xttrvar = dataset(log10(1./Refl_Tot(:,wavelenghtselection )));
Xttrvar.axisscale{2} = 1:length(wavelenghtselection);
Xttrvar.include{2}=Model_par_use{g};
[press,cumpress,rmsecv,rmsec,cvpred,misclassed] = crossval(Xttrvar(1:150,:),Y,'sim',Cow_ID_Tot(1:150),ncomp,Model_par_optionspls{g});
% Calculate absolute residuals
for j = 1:ncomp
    absres(:,j) = abs(cvpred(:,j)-Y);
end
% select number of latent variables based on F test
[cumpressmin,h_star] = min(cumpress);
Hlv = 1;
h = 1;
F_amm = ftest(0.25,length(Y),length(Y));
while Hlv == 1;
    F = cumpress(h)/cumpress(h_star);
    if F < F_amm
        Hlv = 0;
    else
        h = h+1;
    end
end
nsel = h;
clear h Hlv F h_star F_amm;
modelvarpls = pls(Xttrvar(1:150,:),Y,nsel,Model_par_optionspls{g});
validvarpls = pls(Xttrvar(151:end,:),Yt,modelvarpls,Model_par_optionspls{g});
Lactose_Refl2 = cvpred(:,nsel);
Lactose_Refl2t = cell2mat(validvarpls.pred)';
figure(99),subplot(2,3,6),plot(Y,cvpred(:,nsel)','b*',Yt,cell2mat(validvarpls.pred)','go',[0,10],[0,10],'k-')
legend({['RMSE_C_V = ' sprintf('%0.3f',rmsecv(nsel))],['RMSE_P = ' sprintf('%0.3f',sqrt(nanmean((cell2mat(validvarpls.pred)-Yt).^2)))]},'location','NorthWest','FontSize',12);
title(['Reflectance - LV = ' num2str(nsel)],'FontSize',14)
xlabel('Lactose (%w/w) by reference','FontSize',14)
ylabel('Lactose (%w/w) by sensor','FontSize',14)
axis([3.5,5.5,3.5,5.5])
figure(100),subplot(2,3,6),plot([1:length(cvpred(:,nsel))],cvpred(:,nsel)-Y,[1:length(Yt)]+length(cvpred(:,nsel)),(cell2mat(validvarpls.pred)-Yt),[0,1300],[0,0],'k-')
line(I(:)',J(1:length(I(:))),'color','k','LineWidth',2);
ylabel('Error (%w/w)    [Y_{Pred} - Y_{Ref}]','FontSize',14)
title('Reflectance - Lactose','FontSize',14)
xlabel('Samples','FontSize',14)
ylim([-1,1])
xlim([0,1300])

%% Plot spectra

figure(),plot(Trans_Tot','linewidth',2)
set(gca,'XTick',[24:33:256],'XTickLabel',[1000:100:1700]);
xlabel('Wavelength (nm)','FontSize',14)
ylabel('Transmittance (%)','FontSize',14)
title('Transmittance','FontSize',18)
xlim([0,256])

figure(),plot(Refl_Tot'*100,'linewidth',2)
set(gca,'XTick',[24:33:256],'XTickLabel',[1000:100:1700]);
xlabel('Wavelength (nm)','FontSize',14)
ylabel('Reflectance (%)','FontSize',14)
title('Reflectance','FontSize',18)
xlim([0,256])













