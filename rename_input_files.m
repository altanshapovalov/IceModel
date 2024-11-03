% Renaming input files for the new icemodel
% Created Jan 27, 2024; Maxim Altan-Lu Shapovalov

% PURPOSE:rename existing input files of met (AWS data) and albedos (MODIS, MAR,
% RACMO, MERRA-2) used for old icemodel (pre-2024) and out them into
% icemodel>input met and userdata so that icemodel_new can read them and run

% *************
% IMPORTANT: this code (1.1 and 1.2) was run on Jan 27, 2024 and the input files
% were successfully transfered (to icemodel>input>met/userdata). This code 
% should not be run again to prevent overwriting of functional met files; use it
% only for reference.
% *************

clear; clc

%%  0 ) Specify parameters for file extraction
% sitename: PROMICE AWS site
% userdata: from which source will I replace AWS's albedo
sitename = 'kanl';        % options: 'kanm' 'kanl' 'kanu'
userdata = 'mar';       % options: 'modis','racmo','merra','mar'
simyears = 2010:2021;     % study period (simulation years)
forcings = sitename;      % the forcing will always match the sitename

mp = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/'; % mp = max's path



%% 1) Obtain the specific input file used for old Icemodel


% %% 1.1) met
% % a) specify path where all met files are saved
% loc = [mp 'drive_python/outputs/runoff_temp_feedback/' sitename '/read_Promice_AWS/mets_2009_2022.mat'];
% load(loc); % load the met files from that path (loc=location)
% 
% % b) rename and put in new place
% for yr = simyears
%     met = eval(['met' num2str(yr)]); % rename the file of a given year to met
%     outloc = [mp 'input/met/']; % output location (save to this new place)
%     newname = ['met_' sitename '_' forcings '_' num2str(yr) '_1hr.mat'];
%     nnal = [outloc newname]; % new name and location (for met file)
%     save(nnal,'met'); % save in the new location
% end
% 
% % c) check if files are the same
% clear met; % clear any met that could be saved
% testyear = 2012; % specify year we wanna test
% newmet = ['met_' sitename '_' forcings '_' num2str(testyear) '_1hr.mat']; %name
% getfromhere = [outloc newmet]; % extract the new met file for comparison
% load(getfromhere); % load the new met file
% % specify the og file we're comparing it to:
% og_met = eval(['met' num2str(testyear)]);
% % visual comparison:
% plot(og_met.Time,og_met.tair); hold on; pause; plot(met.Time,met.tair)
% % numerical confirmation of match:
% howsimilar = corr(og_met.tair,met.tair); disp(howsimilar);
% % if 1.0, then it's identical = file transfered correctly!


% %% 1.2) albedos (userdata)
% % a) specify path where albedo files are saved
% loc = [mp 'drive_python/outputs/manuscript/albedo/' userdata '/' sitename '/' userdata '_albs_2009_2022.mat'];
% load(loc); % load the met files from that path (loc=location)
% 
% % b) rename and put in new place
% for yr = simyears
%     Data = eval([userdata num2str(yr)]); % rename the file of a given year to met
%     outloc = [mp 'input/userdata/']; % output location (save to this new place)
%     newname = [userdata '_' sitename '_' num2str(yr) '.mat'];
%     nnal = [outloc newname]; % new name and location (for met file)
%     save(nnal,'Data'); % save in the new location
% end
% 
% % c) check if files are the same
% clear Data; % clear any met that could be saved
% testyear = 2012; % specify year we wanna test
% newData = [userdata '_' sitename '_' num2str(testyear) '.mat']; %name
% getfromhere = [outloc newData]; % extract the new met file for comparison
% load(getfromhere); % load the new met file
% % specify the og file we're comparing it to:
% og_Data = eval([userdata num2str(testyear)]);
% % visual comparison:
% % plot(og_Data.Time,og_Data.albedo); hold on; pause; plot(Data.Time,Data.albedo)
% % numerical confirmation of match:
% howsimilar = corr(og_Data.albedo,Data.albedo); 
% disp([sitename ', ' userdata ' albedo correlation: ' num2str(howsimilar)]);
% % if 1.0, then it's identical = file transfered correctly!
% 
% 


