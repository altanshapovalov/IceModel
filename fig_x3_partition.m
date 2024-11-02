% Figure 3 – Average Albedos for all models at all AWS 2010-2021
% by Maxim Altan-Lu Shapovalov for a thesis manuscript
% Nov 2023

% Contents:
% 0) Preamble
% 1) Plot average albedos (together)
% 2) Pre-processing
    % 2.1) Organize all albedos
    % 2.2) Calculate mins, maxs, and averages

% NOTE: Pre-processing has code that prepared the files that get loaded in the
% beginning of section 1. It's best to just run section 1 and get the figure;
% sections 2.1 and 2.2 are more for seeing how those files were created

%% Preamble


clear; clc; % clear workspace and command window for neatness
% make future paths shorter with a shortcut (mp = max's path)
mp = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/'; 
data = {'modis' 'merra' 'mar' 'racmo'}; % albedo sources

% modidx = 1;


% % specify AWSs
% aws = {'kanl' 'kanm' 'kanu'};
aws = {'kanl'};
% aws = {'kanm'};
% aws = {'kanl','kanm'};
% aws = {'kpcl' 'kpcu'};
% aws = {'nukl' 'nuku'};
% aws = {'qasl' 'qasu'};
% aws = {'scol' 'scou'};
% aws = {'thul' 'thuu'};
% aws = {'upel' 'upeu'};
% aws = {'kanl','kanm','kanu','kpcl','kpcu','nukl','nuku','qasl','qasu','scol',...
%     'scou', 'thul', 'thuu', 'upel', 'upeu'};

% establish colors for the plotting
colors = ({[0 0 0]; % modis (black)
    [0.9 0.5 0.1]; % merra (red-orange)
    [0 0.4470 0.7410]; % mar (dark blue)
%     [0.4010 0.8450 0.9330]})'; % og racmo (light blue)
    [0.6350 0.0780 0.1840]})'; % new racmo (maroon)

% Establish albedo thresholds:
glacier_ice_thresh = 0.55; % Ryan et al., 2019
glacier_ice_thresh_merra = 0.60; % Cullather et al., 2014

% save figs here
% FolderPath = [mp 'drive_python/figures/29. new aws/25. partition tables/ice/figs/vs elevation/'];

% logging output info (attencion: concatenates!)
% DiaryPath = [mp 'drive_python/figures/29. new aws/25. partition tables/'];
% DiaryName = 'all_aws_partition_tables.txt';
% diary(fullfile(DiaryPath, DiaryName));


%% Partitioning

% Define the number of rows in the table
num_rows = length(data); % each row for each model minus modis
colN = nan(num_rows, 1); % Column with double
% prep individual arrays for creation of all_aws_partition table (next section):
col1 = nan(num_rows,length(aws));
col2 = col1; col3 = col1; col4 = col1; col5 = col1; col6 = col1; col7 = col1;
col8 = col1; col9 = col1; col10 = col1; col11 = col1; col12 = col1;
col13 = col1; col14 = col1; col15 = col1;



for a = 1:length(aws) % 1st, loop thru weather stations
    met = aws{a}; % choose current aws
    alb(8760,length(data)) = NaN; alb(:,:) = NaN; % prep for storage
    % create empty table for storage of partition data
    data_table = table(colN,colN,colN,colN,colN,...
        colN,colN,colN,colN,colN,...
        colN,colN,colN,colN,colN,... % prep 15 NaN columns
        'VariableNames',{'P1_dA','P2_dA','P3_dA','P4_dA','P5_dA',... % alb differences
        'P1_dM','P2_dM','P3_dM','P4_dM','P5_dM',... % melt differences
        'P1','P2','P3','P4','P5'},...% period lengths
        'RowNames',{'MODIS','MERRA','MAR','RACMO'});

    % Get albedo data
    for b = 1:length(data) % 2nd, loop thru diff models
        % i) load in file (will be named "avg")
        if b == 1 % if modis, use new albedo (mcd43a43)
            load([mp 'drive_python/outputs/manuscript/albedo/new_aws/mcd43a3/single pixel/' met '_' data{b} '.mat'])
        else % for other data sources
            load([mp 'drive_python/outputs/manuscript/albedo/new_aws/' met '_' data{b} '.mat'])
        end
        % ii) put data in a structure
        alb(:,b) = avg.albedo;
    end

    % Get melt data
    for b = 1:length(data) % 2nd, loop thru diff models
        % i) load in file (will be named "avg")
        if b==1 % if modis, use mcd43a3
            load([mp 'drive_python/outputs/manuscript/melt/new_aws/mcd43a3/single pixel/' met '_' data{b} '.mat'])
        else
            load([mp 'drive_python/outputs/manuscript/melt/new_aws/' met '_' data{b} '.mat'])
        end
        % ii) save individually
        m50(:,b) = avg.melt50;
    end

    % rename data for code consistency (i.e., chunks of code below were written
    % with vars avgsdata, avg50, and tiempo1 in mind)
    avgdata = alb; avg50 = m50; tiempo1 = avg.Time;

    % Begin partitioning
    % i) find when modis alb shifts snow to ice and vice versa
    for cc = 1:length(data)

        % Give merra-2 the threshold that it uses (0.6 vs 0.55, as in literature)
        if cc == 2 % if merra
            idxs = find(avgdata(:,cc) <= glacier_ice_thresh_merra); % all instances when modis-alb is 0.6 or lower
        else % for all others use regular 0.55 threshold for glacier ice
            idxs = find(avgdata(:,cc) <= glacier_ice_thresh); % all instances when modis-alb is 0.55 or lower
        end

        if ~isempty(idxs) % if there is at least 1 instance when EGI (exposed glacial ice)
            % a) transition days (snow, snowline, ice)
            idx1 = min(idxs); % id 1st instance
            idx2 = max(idxs); % and last instance
            % save when MODIS specifically has transitions because will compare
            % other models to it:
            if cc == 1 % if MODIS
                modis_idx1 = idx1; modis_idx2 = idx2; % will compare to these
            end

            % b) calculate length of periods
            % Concept (when comparing modis to another model):
            % t1 – row of the 1st model that exposes glacier ice (EGI)
            % t2 - row of the 2nd model EGI
            % t3 – row of the 1st model's last EGI
            % t4 – row of the 2nd model's last EGI
            % t0 – first row of the year
            % t_end – last row/hour of the year
            % So, doesn't matter which model EGI 1st, this will track it
            % regardless of model name.

            if cc==1 % if modis
                % actually, don't do anything
%                 mod_alb_sno1  = median(avgdata(1:modis_idx1,cc)); % avg modis snow albedo pre EGI
%                 mod_alb_ice   = median(avgdata(modis_idx1:modis_idx2,cc)); % avg modis ice albedo
%                 mod_alb_sno2  = median(avgdata(modis_idx2:length(avgdata),cc)); % avg modis snow albedo post EGI
%                 mod_melt_sno1 = avg50(modis_idx1,cc); % modis melt at 1st EGI
%                 mod_melt_ice  = avg50(modis_idx2,cc); % modis melt at last EGI
%                 mod_melt_sno2 = avg50(length(avgdata),cc); % modis melt at last hour of the year
            else % for other models
                % c) establish the transitions (t's)
                % t1 & t2:
                if modis_idx1 < idx1 % if modis is 1st to EGI
                    t1 = modis_idx1; t2 = idx1;
%                     disp(['MODIS EGI before ' data{cc}])
                elseif modis_idx1 > idx1 % if modis is 2nd to EGI
                    t1 = idx1; t2 = modis_idx1;
%                     disp([data{cc} ' EGI before MODIS'])
                end
                % t3 & t4:
                if modis_idx2 < idx2 % if modis is 1st to cover EGI
                    t3 = modis_idx2; t4 = idx2;
%                     disp(['MODIS covers GI before ' data{cc}])
                elseif modis_idx2 > idx2 % if modis is 2nd to cover EGI
                    t3 = idx2; t4 = modis_idx2;
%                     disp([data{cc} ' covers GI before MODIS'])
                end
                % Attempting to know the EGI of each model:
%                 if cc==1 % if modis
%                     disp(['EGI (in days) of MODIS: ' modis_idx2 - modis_idx1])
%                 else
%                     disp(['EGI (in days) of MODIS: ' modis_idx2 - modis_idx1])
%                 end
                t0 = 1; % first row
                t_end = length(avgdata); % last row
    
                % d) establish the periods (P's) (no bias)
%                 P1 = t1 - 1; % both snow
%                 P2 = t2 - t1; % 1 model is snow, other is already EGI
%                 P3 = t3 - t2; % both models are EGI
%                 P4 = t4 - t3; % 1 still EGI, other already snow
%                 P5 = t_end - t4; % both snow once again
%                 disp(['By how much for P2: ' num2str(P2/24)])
%                 disp(['By how much for P3: ' num2str(P3/24)])
%                 disp(['By how much for P4: ' num2str(P4/24)])
%                 disp(' ')
                % d) establish the periods (P's) (WITH bias)
%                 P1 = idx1 - 1; % how much snow a model shows
                P1 = t1 - 1; % when both modis and model are snow
                P2 = idx1 - modis_idx1; % model minus modis. If positive, that 
                % means model's idx (i.e., row, i.e., day) of EGI is later than 
                % MODIS, meaning it exposes EGI later
%                 P3 = idx2 - idx1; % total EGI for a given model
                P3 = t3 - t2; % both models are EGI
                P4 = idx2 - modis_idx2; % if positive, model covers EGI later than modis
%                 P5 = t_end - idx2; % again, each model's snow representation
                P5 = t_end - t4; % both snow once again

    
                % e) get avg albedo for each period
                p1_a = median(avgdata(t0:t1,cc));
                p2_a = median(avgdata(t1:t2,cc));
                p3_a = median(avgdata(t2:t3,cc));
                p4_a = median(avgdata(t3:t4,cc));
                p5_a = median(avgdata(t4:t_end,cc));
                % get it for modis specifically for comparison
                p1_mod_a = median(avgdata(t0:t1,1));
                p2_mod_a = median(avgdata(t1:t2,1));
                p3_mod_a = median(avgdata(t2:t3,1));
                p4_mod_a = median(avgdata(t3:t4,1));
                p5_mod_a = median(avgdata(t4:t_end,1));
                % save all the albs
                albz(1,a) = p3_mod_a; % modis in 1st row
                if cc>1 % for rows 2-4, other models
                    albz(cc,a) = p3_a;
                end

                % f) get melts for each period
                p1_m = avg50(t1,cc);
                p2_m = avg50(t2,cc);
                p3_m = avg50(t3,cc);
%                 disp([met ' ' data{cc} ' melt: ' num2str(p3_m)])
                p4_m = avg50(t4,cc);
                p5_m = avg50(t_end,cc);
                % for modis:
                p1_mod_m = avg50(t1,1);
                p2_mod_m = avg50(t2,1);
                p3_mod_m = avg50(t3,1);
                p4_mod_m = avg50(t4,1);
                p5_mod_m = avg50(t_end,1);
%                 disp([met ' ' data{1} ' melt: ' num2str(p3_mod_m)])
                % save all the melts
                meltz(1,a) = p5_mod_m; % modis in 1st row
                if cc>1 % for rows 2-4, other models
                    meltz(cc,a) = p5_m;
                end

                % g) Calculate albedo and melt differences (regular and RMSEs):
                p1_dA = p1_a - p1_mod_a; % albedos
                p2_dA = p2_a - p2_mod_a;
                p3_dA = p3_a - p3_mod_a;
                p4_dA = p4_a - p4_mod_a;
                p5_dA = p5_a - p5_mod_a;
                % THIS^ gives actual albedo differences (not as percentages)
                % calculate cumulative melt differences:
                p1_dM = p1_m - p1_mod_m; % melts
                p2_dM = p2_m - p2_mod_m;
                p3_dM = p3_m - p3_mod_m;
                p4_dM = p4_m - p4_mod_m;
                p5_dM = p5_m - p5_mod_m;
                % check each EGI melt calculation:
%                 disp([met ' ' data{cc}])
%                 disp(['p3_m : ' num2str(p3_m)])
%                 disp(['p3_mod_m : ' num2str(p3_mod_m)])
%                 disp(['p3_dM : ' num2str(p3_dM)])
%                 disp(' ')
                % calculate melt made in individual period:
                p1_dM = p1_m - p1_mod_m;
                p2_dM = (p2_m - p1_m) - (p2_mod_m - p1_mod_m);
                p3_dM = (p3_m - p2_m) - (p3_mod_m - p2_mod_m);
                p4_dM = (p4_m - p3_m) - (p4_mod_m - p3_mod_m);
                p5_dM = (p5_m - p4_m) - (p5_mod_m - p4_mod_m);


                % AS PERCENTAGES:
                % g) Calculate albedo and melt differences (regular and RMSEs):
%                 p1_dA = ((p1_a - p1_mod_a)/p1_mod_a)*100; % albedos
%                 p2_dA = ((p2_a - p2_mod_a)/p2_mod_a)*100;
%                 p3_dA = ((p3_a - p3_mod_a)/p3_mod_a)*100;
%                 p4_dA = ((p4_a - p4_mod_a)/p4_mod_a)*100;
%                 p5_dA = ((p5_a - p5_mod_a)/p5_mod_a)*100;
                % THIS RIGHT HERE (^) gives albedo differences as perecentages
                % cumulative melts:
%                 p1_dM = ((p1_m - p1_mod_m)/p1_mod_m)*100; 
%                 p2_dM = ((p2_m - p2_mod_m)/p2_mod_m)*100;
%                 p3_dM = ((p3_m - p3_mod_m)/p3_mod_m)*100;
%                 p4_dM = ((p4_m - p4_mod_m)/p4_mod_m)*100;
%                 p5_dM = ((p5_m - p5_mod_m)/p5_mod_m)*100;
                % individual melts as perecntages:(nein)
%                 p3_dM = p3_dM/p3_mod_m * 100;



                % i) put into a table unqiue to the model
                data_table{cc,1} = p1_dA; % dAlb
                data_table{cc,2} = p2_dA;
                data_table{cc,3} = p3_dA;
                data_table{cc,4} = p4_dA;
                data_table{cc,5} = p5_dA;
                data_table{cc,6} = p1_dM; % dMelt
                data_table{cc,7} = p2_dM;
                data_table{cc,8} = p3_dM;
                data_table{cc,9} = p4_dM;
                data_table{cc,10} = p5_dM;
                data_table{cc,11} = P1; % periods
                data_table{cc,12} = P2;
                data_table{cc,13} = P3;
                data_table{cc,14} = P4;
                data_table{cc,15} = P5;

                % j) allocate to unqiue vars for the 'mega' table
                col1(cc,a) = p1_dA;
                col2(cc,a) = p2_dA;
                col3(cc, a) = p3_dA;
                col4(cc, a) = p4_dA;
                col5(cc, a) = p5_dA;
                col1(1,a) = p1_mod_a;
                col2(1,a) = p2_mod_a;
                col3(1,a) = p3_mod_a;
                col4(1,a) = p4_mod_a;
                col5(1,a) = p5_mod_a;
                col6(cc, a) = p1_dM;
                col7(cc, a) = p2_dM;
                col8(cc, a) = p3_dM;
                col9(cc, a) = p4_dM;
                col10(cc, a) = p5_dM;
                col6(1,a) = p1_mod_m;
                col7(1,a) = p2_mod_m;
                col8(1,a) = p3_mod_m;
                col9(1,a) = p4_mod_m;
                col10(1,a) = p5_mod_m;
                col11(cc, a) = P1; 
                col12(cc, a) = P2;
                col13(cc, a) = P3;
                col14(cc, a) = P4;
                col15(cc, a) = P5;

                % k) separately save the avg (median) melts (M), albs (A), and total EGIs (P)
%                 AvgIceInfo_P(cc,a) = P3 / 24; % convert to days
%                 AvgIceInfo_A(cc,a) = 
               
            end % if modis, this conditional above is skipped

        end % if no EGI, loop above is skipped

    end % finished running for all 4 models

    % rename data_table to have a name unique to the AWS
    dt_nom = [met '_partition_data']; % the name itself
    assignin('base',dt_nom,data_table) % create the var in workspace

    clear data_table; % so nothing gets overwritten for different AWSs

end

% Now we have a table for each AWS! Onward to create a 'mega' table that holds
% median/average values of all AWS tables


% %% Alb extraction
% albz(albz==0) = NaN; % get rid of zeros that can skew my results
% aa(1,1) = median(albz(1,:),'omitnan');
% aa(2,1) = median(albz(2,:),'omitnan');
% aa(3,1) = median(albz(3,:),'omitnan');
% aa(4,1) = median(albz(4,:),'omitnan');
% 
% egi_albz = table(aa,'VariableNames', {'Albedo'},...
%     'RowNames', {'MODIS', 'MERRA', 'MAR', 'RACMO'});
% disp('Avg albedo for EGI')
% disp(egi_albz)
% 
% %% Melt Extraction (getting the avg EGI value)
% meltz(meltz==0) = NaN; % get rid of zeros that can skew my results
% mm(1,1) = median(meltz(1,:),'omitnan');
% mm(2,1) = median(meltz(2,:),'omitnan');
% mm(3,1) = median(meltz(3,:),'omitnan');
% mm(4,1) = median(meltz(4,:),'omitnan');
% 
% 
% egi_meltz = table(mm,...
%     'VariableNames', {'Melt'},...
%     'RowNames', {'MODIS', 'MERRA', 'MAR', 'RACMO'});
% disp('Avg melt for EGI')
% disp(egi_meltz)

% %% 
% % a) id columns when RACMO is NaN
% nan_columns = any(isnan(meltz(4,:)), 1);
% 
% % Extract numbers that are NOT in the columns with NaNs
% newmod = meltz(1, ~nan_columns);
% newrac = meltz(4, ~nan_columns);
% newdiff = newrac - newmod;
% racmed = median(newrac) 
% modmed = median(newmod) 
% diffmed = mean(newdiff)





% %% Мега таблица
% prep arrays for preallocation (increases speed)
c1 = nan(num_rows,1); 
c2 = c1; c3 = c1; c4 = c1; c5 = c1; c6 = c1; c7 = c1; c8 = c1; c9 = c1; 
c10 = c1; c11 = c1; c12 = c1; c13 = c1; c14 = c1; c15 = c1;
% get avg values for individual arrays
for dd = 1:length(data)
    c1(dd, 1) = median(col1(dd, :), 'omitnan');
    c2(dd, 1) = median(col2(dd, :), 'omitnan');
    c3(dd, 1) = median(col3(dd, :), 'omitnan');
    c4(dd, 1) = median(col4(dd, :), 'omitnan');
    c5(dd, 1) = median(col5(dd, :), 'omitnan');
    c6(dd, 1) = median(col6(dd, :), 'omitnan');
    c7(dd, 1) = median(col7(dd, :), 'omitnan');
    c8(dd, 1) = median(col8(dd, :), 'omitnan');
    c9(dd, 1) = median(col9(dd, :), 'omitnan');
    c10(dd, 1) = median(col10(dd, :), 'omitnan');
    c11(dd, 1) = median(col11(dd, :), 'omitnan');
    c12(dd, 1) = median(col12(dd, :), 'omitnan');
    c13(dd, 1) = median(col13(dd, :), 'omitnan');
    c14(dd, 1) = median(col14(dd, :), 'omitnan');
    c15(dd, 1) = median(col15(dd, :), 'omitnan');
end

all_aws_partition = table(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,...
    'VariableNames', {'P1_dA','P2_dA','P3_dA','P4_dA','P5_dA',... % alb differences
    'P1_dM','P2_dM','P3_dM','P4_dM','P5_dM',... % melt differences
    'P1','P2','P3','P4','P5'},...% period lengths
    'RowNames', {'MODIS', 'MERRA', 'MAR', 'RACMO'});
disp('Partition table for albedo and melt differences for all AWS:')
disp(all_aws_partition)




%% RMSE tables

% ALBEDO:
% prepare for storage
ra1 = nan(num_rows,1); ra2 = ra1; ra3 = ra1; ra4 = ra1; ra5 = ra1; ra6 = ra1;
ra1p = ra1; ra2p = ra1; ra3p = ra1; ra4p = ra1; ra5p = ra1; % p is for "%"
% RMSEs for each period
for gg = 1:length(data)
    ra1(gg,1) = sqrt(mean(c1(gg,1).^2,'omitnan'));
    ra2(gg,1) = sqrt(mean(c2(gg,1).^2,'omitnan'));
    ra3(gg,1) = sqrt(mean(c3(gg,1).^2,'omitnan'));
    ra4(gg,1) = sqrt(mean(c4(gg,1).^2,'omitnan'));
    ra5(gg,1) = sqrt(mean(c5(gg,1).^2,'omitnan'));
end
% Total RMSE for each model
for gg = 1:length(data) 
    ra6(gg,1) = ra1(gg,1) + ra2(gg,1) + ra3(gg,1) + ra4(gg,1) + ra5(gg,1);
end
% RMSEs as percentages (to really help partition)
for gg = 1:length(data)
    ra1p(gg,1) = (ra1(gg,1) / ra1(1,1)) * 100; % convert to percent
    ra2p(gg,1) = (ra2(gg,1) / ra2(1,1)) * 100;
    ra3p(gg,1) = (ra3(gg,1) / ra3(1,1)) * 100;
    ra4p(gg,1) = (ra4(gg,1) / ra4(1,1)) * 100;
    ra5p(gg,1) = (ra5(gg,1) / ra5(1,1)) * 100;
end
% create the table with the values
alb_rmse_table = table(ra1,ra1p,ra2,ra2p,ra3,ra3p,ra4,ra4p,ra5,ra5p,ra6,...
    'VariableNames',{'P1','P1(%)','P2','P2(%)','P3','P3(%)',...
    'P4','P4(%)','P5','P5(%)','Total'},...% period lengths
    'RowNames',{'MODIS', 'MERRA', 'MAR', 'RACMO'});
disp('RMSE of avg albedo differences for all AWS')
disp(alb_rmse_table)



% MELT:
% prepare for storage
rm1 = nan(num_rows,1); rm2 = rm1; rm3 = rm1; rm4 = rm1; rm5 = rm1; rm6 = rm1;
rm1p = rm1; rm2p = rm1; rm3p = rm1; rm4p = rm1; rm5p = rm1; % p is for "%"
for gg = 1:length(data)
    rm1(gg,1) = sqrt(mean(c6(gg,1).^2,'omitnan'));
    rm2(gg,1) = sqrt(mean(c7(gg,1).^2,'omitnan'));
    rm3(gg,1) = sqrt(mean(c8(gg,1).^2,'omitnan'));
    rm4(gg,1) = sqrt(mean(c9(gg,1).^2,'omitnan'));
    rm5(gg,1) = sqrt(mean(c10(gg,1).^2,'omitnan'));
end
% Total RMSE for each model
for gg = 1:length(data) 
    rm6(gg,1) = rm1(gg,1) + rm2(gg,1) + rm3(gg,1) + rm4(gg,1) + rm5(gg,1);
end
% RMSEs as percentages (to really help partition)
for gg = 1:length(data)
    rm1p(gg,1) = (rm1(gg,1) / rm1(1,1)) * 100; % convert to percent
    rm2p(gg,1) = (rm2(gg,1) / rm2(1,1)) * 100;
    rm3p(gg,1) = (rm3(gg,1) / rm3(1,1)) * 100;
    rm4p(gg,1) = (rm4(gg,1) / rm4(1,1)) * 100;
    rm5p(gg,1) = (rm5(gg,1) / rm5(1,1)) * 100;
end

% create the table with the values
melt_rmse_table = table(rm1,rm1p,rm2,rm2p,rm3,rm3p,rm4,rm4p,rm5,rm5p,rm6,...
    'VariableNames',{'P1','P1(%)','P2','P2(%)','P3','P3(%)',...
    'P4','P4(%)','P5','P5(%)','Total'},...% period lengths
    'RowNames',{'MODIS', 'MERRA', 'MAR', 'RACMO'});
disp('RMSE of avg melt differences for all AWS')
disp(melt_rmse_table)


%% RMSE Albedo bar graph

% ALBEDO:
% prepare for storage
ra1 = nan(num_rows,1); ra2 = ra1; ra3 = ra1; ra4 = ra1; ra5 = ra1; ra6 = ra1;
ra1p = ra1; ra2p = ra1; ra3p = ra1; ra4p = ra1; ra5p = ra1; % p is for "%"
% RMSEs for each period
for gg = 1:length(data)
    ra1(gg,1) = sqrt(mean(c1(gg,1).^2,'omitnan'));
    ra2(gg,1) = sqrt(mean(c2(gg,1).^2,'omitnan'));
    ra3(gg,1) = sqrt(mean(c3(gg,1).^2,'omitnan'));
    ra4(gg,1) = sqrt(mean(c4(gg,1).^2,'omitnan'));
    ra5(gg,1) = sqrt(mean(c5(gg,1).^2,'omitnan'));
end
% Total RMSE for each model
for gg = 1:length(data) 
    ra6(gg,1) = ra1(gg,1) + ra2(gg,1) + ra3(gg,1) + ra4(gg,1) + ra5(gg,1);
end
% RMSEs as percentages (to really help partition)
for gg = 1:length(data)
    ra1p(gg,1) = (ra1(gg,1) / ra6(gg,1)) * 100; % convert to percent
    ra2p(gg,1) = (ra2(gg,1) / ra6(gg,1)) * 100;
    ra3p(gg,1) = (ra3(gg,1) / ra6(gg,1)) * 100;
    ra4p(gg,1) = (ra4(gg,1) / ra6(gg,1)) * 100;
    ra5p(gg,1) = (ra5(gg,1) / ra6(gg,1)) * 100;
end

% MELT:
% prepare for storage
rm1 = nan(num_rows,1); rm2 = rm1; rm3 = rm1; rm4 = rm1; rm5 = rm1; rm6 = rm1;
rm1p = rm1; rm2p = rm1; rm3p = rm1; rm4p = rm1; rm5p = rm1; % p is for "%"
for gg = 1:length(data)
    rm1(gg,1) = sqrt(mean(c6(gg,1).^2,'omitnan'));
    rm2(gg,1) = sqrt(mean(c7(gg,1).^2,'omitnan'));
    rm3(gg,1) = sqrt(mean(c8(gg,1).^2,'omitnan'));
    rm4(gg,1) = sqrt(mean(c9(gg,1).^2,'omitnan'));
    rm5(gg,1) = sqrt(mean(c10(gg,1).^2,'omitnan'));
end
% Total RMSE for each model
for gg = 1:length(data) 
    rm6(gg,1) = rm1(gg,1) + rm2(gg,1) + rm3(gg,1) + rm4(gg,1) + rm5(gg,1);
end
% RMSEs as percentages (to really help partition)
for gg = 1:length(data)
    rm1p(gg,1) = (rm1(gg,1) / rm6(gg,1)) * 100; % convert to percent
    rm2p(gg,1) = (rm2(gg,1) / rm6(gg,1)) * 100;
    rm3p(gg,1) = (rm3(gg,1) / rm6(gg,1)) * 100;
    rm4p(gg,1) = (rm4(gg,1) / rm6(gg,1)) * 100;
    rm5p(gg,1) = (rm5(gg,1) / rm6(gg,1)) * 100;
end

% Reassign var name for easier switch
% baropt = 'alb';
baropt = 'melt';
switch baropt
    case 'alb'
        y1 = ra1;
        y2 = ra2;
        y3 = ra3;
        y4 = ra4;
        y5 = ra5;
        rect_height = 0.16; % top bound for rectangle
        y_ticks = 0:0.02:rect_height;
    case 'melt'
        y1 = rm1;
        y2 = rm2;
        y3 = rm3;
        y4 = rm4;
        y5 = rm5;
        rect_height = 0.7; % top bound for rectangle
        y_ticks = 0:0.1:rect_height;
end

% Define the names of the groups
group_names = {'P1', 'P2', 'P3', 'P4', 'P5'};
% Define bar colors for each group
colors = {[0.9 0.5 0.1],... % merra2 (orange)
    [0 0.4470 0.7410],... % mar (dark blue)
    [0.6350 0.0780 0.1840]}; % racmo (dark red)
% control bar width (1 makes them perfectly flush)
bw = 1;
% Define background colors for each group
background_colors = [
    0.9 0.9 0.9;  % Light gray for groups P1 and P5
    0.8 0.8 0.8;  % Slightly darker gray for groups P2 and P4
    0.7 0.7 0.7;  % darkest for P3
    0.8 0.8 0.8;  % Darker gray for rectangle 4
    0.9 0.9 0.9]; % Lightest gray for rectangle 5];  % Darkest gray for group P3
% Define the x-axis positions for the rectangles
x_positions = [0 4 8 12 16 20];

figure
% Plot proxy objects for legend with specified colors
for i = 1:numel(colors)
    bar(NaN, 'FaceColor',colors{i}); % Invisible bar for legend
    hold on;
end
% Create background areas between groups
% Iterate over each rectangle
for i = 1:size(background_colors, 1)
    % Define the x and y coordinates of the rectangle
    x_rect = [x_positions(i) x_positions(i+1) x_positions(i+1) x_positions(i)];
    y_rect = [0 0 rect_height rect_height]; % Rectangle extends from y=0 to y=1
    
    % Create the rectangle patch
    patch(x_rect, y_rect, background_colors(i, :), 'EdgeColor', 'none');
    hold on;
end

% doing it twice just so that the legend can show the model colors properly
for i = 1:numel(group_names) % for each period
    % Extract data for the current group
    yy = eval(['y' num2str(i)]); % get current data
    y = yy(:,1); % specific period
    y = y(2:4,1); % ignore the MODIS 1st row NaN
    % Define x-axis positions for the bars in the current group
    x_positions = (1:3) + (i - 1) * 4; % Adjusted position for each group
    % Create bars for the current group with specified color
    for j = 1:3
        bar(x_positions(j),y(j),bw,'FaceColor',colors{j}); 
        hold on;
    end
end

% Set x-axis tick positions and labels
xticks(2:4:20);
xticklabels(group_names);

xlim([0 20])

% Add grid lines since "grid on" doesn't work with rectangles
% Define the x-axis and y-axis tick positions
% x_ticks = 0:2:20;
% Get the x-axis and y-axis limits
x_limits = xlim;
y_limits = ylim;
% line color (for grid)
lc = [0.5 0.5 0.5];
% Plot horizontal lines at each y tick mark
for y_tick = y_ticks
    line(x_limits, [y_tick y_tick], 'Color', lc, 'LineStyle', '-', 'LineWidth', 0.1,'HandleVisibility', 'off');
end

% Add labels and title
switch baropt
    case 'alb'
        ylabel('Albedo [-]');
%         title("RMSE of each model's albedo relative to MODIS albedo");
%         legend('MERRA-2', 'MAR', 'RACMO','Snow','Snowline','Ice','Position',[0.158, 0.777, 0.1, 0.11])
    case 'melt'
        ylabel('Melt [m]');
%         title("RMSE of melt for each model's albedo relative to MODIS");
%         legend('MERRA-2', 'MAR', 'RACMO','Snow','Snowline','Ice','Position',[0.158, 0.777, 0.1, 0.11])
end
xlabel('Periods');
grid on;


%% Box plots for albedo and melt biases

% Reassign var name for easier switch
% baropt = 'alb';
% baropt = 'melt';
baropt = 'period';
switch baropt
    case 'alb'
        y1 = col1;
        y2 = col2;
        y3 = col3;
        y4 = col4;
        y5 = col5;
        rh_t = 0.2; % top bound for rectangle
        rh_b = -0.2; % bottom bound for rectangle
        y_ticks = rh_b:0.05:rh_t;
    case 'melt'
        y1 = col6;
        y2 = col7;
        y3 = col8;
        y4 = col9;
        y5 = col10;
        rh_t = 0.7; % top bound for rectangle
        rh_b = -0.7; % bottom bound for rectangle
        y_ticks = rh_b:0.1:rh_t;
    case 'period'
        y1 = col11/24;
        y2 = col12/24;
        y3 = col13/24;
        y4 = col14/24;
        y5 = col15/24;
%         rect_height = 0.6; % top bound for rectangle
        rh_t = 200; % top bound for rectangle
        rh_b = -20; % bottom bound for rectangle
        y_ticks = rh_b:20:rh_t;

end

% Define the names of the groups
group_names = {'P1', 'P2', 'P3', 'P4', 'P5'};
% Define bar colors for each group
colors = {[0.9 0.5 0.1],... % merra2 (orange)
    [0 0.4470 0.7410],... % mar (dark blue)
    [0.6350 0.0780 0.1840]}; % racmo (dark red)
% control box width (1 makes them perfectly flush)
bw = 0.8;
% Define background colors for each group
background_colors = [
    0.9 0.9 0.9;  % Light gray for groups P1 and P5
    0.8 0.8 0.8;  % Slightly darker gray for groups P2 and P4
    0.7 0.7 0.7;  % darkest for P3
    0.8 0.8 0.8;  % Darker gray for rectangle 4
    0.9 0.9 0.9]; % Lightest gray for rectangle 5];  % Darkest gray for group P3
% Define the x-axis positions for the rectangles
x_positions = [0 4 8 12 16 20];

fig = figure;
% Plot proxy objects for legend with specified colors
for i = 1:numel(colors)
    bar(NaN, 'FaceColor',colors{i}); % Invisible bar for legend
    hold on;
end
% Create background areas between groups
% Iterate over each rectangle
for i = 1:size(background_colors, 1)
    % Define the x and y coordinates of the rectangle
    x_rect = [x_positions(i) x_positions(i+1) x_positions(i+1) x_positions(i)];
    y_rect = [rh_b rh_b rh_t rh_t]; % Rectangle extends from y=0 to y=1
    
    % Create the rectangle patch
    patch(x_rect, y_rect, background_colors(i, :), 'EdgeColor', 'none');
    hold on;
end

% Add grid lines since "grid on" doesn't work with rectangles
% Define the x-axis and y-axis tick positions
x_ticks = 0:2:20;
% Get the x-axis and y-axis limits
x_limits = xlim;
y_limits = ylim;
% line color (for grid)
lc = [0.5 0.5 0.5];
% Plot horizontal lines at each y tick mark
for y_tick = y_ticks
    line(x_limits, [y_tick y_tick], 'Color', lc, 'LineStyle', '-', 'LineWidth', 0.1,'HandleVisibility', 'off');
end
% mark the 0 in a more special format
line(x_limits, [0 0], 'Color', [0.15 0.15 0.15], 'LineStyle', '-', 'LineWidth', 0.65);

for i = 1:numel(group_names) % for each period
    % Extract data for the current group
    yy = eval(['y' num2str(i)]); % get current data
    y = yy(:,:); % specific period
    y = y(2:4,:); % ignore the MODIS 1st row NaN
    % Define x-axis positions for the bars in the current group
    x_positions = (1:3) + (i - 1) * 4; % Adjusted position for each group
    % Create boxplots for the current group with specified color
    for j = 1:3
        h = boxplot(y(j,:),'Positions',x_positions(j),'Colors',colors{j},...
            'Widths',bw,'BoxStyle','outline','Symbol','.','OutlierSize',8);
        boxes = fig.Children.Children(1,1).Children(1:7); % not sure, but get handle of box parts?
        for jj = 1:length(boxes) % draw a colored patch behind each box
            patch(boxes(jj).XData,boxes(jj).YData,colors{j},'FaceAlpha',.5,'EdgeAlpha',0);
        end

        hold on;
        % Extract handles to the individual box graphics objects
        box_handles = findobj(h, 'Tag', 'Box');
        box_handles = findobj(h, 'Tag', 'Box');
        outlier_handles = findobj(h, 'Tag', 'Outliers');
        model_color_index = rem(j - 1, numel(colors)) + 1;
        % Extract handles to the outliers for the current box
        outliers_in_box = outlier_handles(1);
        % Loop over outliers within the current box and assign color
        for k = 1:numel(outliers_in_box)
            outliers_in_box(k).MarkerEdgeColor = colors{model_color_index};
        end
    end
end

% Set x-axis tick positions and labels
xticks(2:4:20);
xticklabels(group_names);

xlim([0 20])

% Add labels and title
switch baropt
    case 'alb'
        ylabel('Albedo [-]');
%         title("Bias of each model's albedo relative to MODIS albedo");
        legend('MERRA-2', 'MAR', 'RACMO','Snow','Snowline','Ice','Position',[0.158, 0.151, 0.1, 0.11])
    case 'melt'
        ylabel('Melt [m]');
%         title("Bias of melt for each model's albedo relative to MODIS");
        legend('MERRA-2', 'MAR', 'RACMO','Snow','Snowline','Ice','Position',[0.158, 0.151, 0.1, 0.11])
    case 'period'
        ylabel('Days');
%         title("Average duration of each period");
        legend('MERRA-2', 'MAR', 'RACMO','Snow','Snowline','Ice','Position',[0.158, 0.151, 0.1, 0.11])
end
xlabel('Periods');
ylim([rh_b rh_t]) % adjust y axis limit based on backrgound rectangle height
grid on; 













