% read MERRA-2 files (M2T3NXGLC) for albedo (as well as surface temp and
% incoming shortwave radiation.

% created on aug 17, 2023

%%
clear; clc
mp = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/'; % mp = max's path

% CHANGE THIS FOR DIFFERENT WEATHER STATION
% aws = 'KANU'; % options: KANM, KANL, KANU

% aws = 'kanl';
aws = 'kanm';
% aws = 'kpcl'; % 1
% aws = 'kpcu'; % 2
% aws = 'nukl'; % 3
% aws = 'nuku'; % 4
% aws = 'qasl'; % 5
% aws = 'qasu'; % 6
% aws = 'scol'; % 7
% aws = 'scou'; % 8
% aws = 'thul'; % 9
% aws = 'thuu'; % 10
% aws = 'upel'; % 11
% aws = 'upeu'; % 12

% % establish the start year (sy)
% if strcmp(aws,'KANM') || strcmp(aws,'KANL') % KANL and KANM start with 2009
%     sy = 2009;
% elseif strcmp(aws,'KANU') % KANU starts with 2010
%     sy = 2009; %changed from 2010 to 2009
% else
%     error([aws ' is NOT a valid weather station on the K-transect (options: KANL, KANM, KANU)']);
% end


%% 1) get files ===========================================================
% Read the names of nc-files in the folder 
pth = [mp 'drive_python/input/MERRA2/glc_select/' aws '/*.nc'];
ncfiles      = dir(pth);   % get basic info
Nfiles       = length(ncfiles);  % Obtain total number of .nc files 

file_name    = ncfiles(1).name;  % extract 1st filename
sy           = str2double(ncfiles(1).name(28:31)); % first available year (start year)
ey           = str2double(ncfiles(end).name(28:31)); % end year
files_time   = length(ncread(file_name,'time')); % if 24 - hourly
ncdisp(file_name); % to see ALL info of .nc file



%% 2) extract vars and put in mr2 structure ===============================

% add time and convert
startDateTime = datetime((sy),1,1,0,0,0); % jan 1 of 1st year
endDateTime = datetime((ey), 12, 31, 23, 0, 0); % dec 31
mr2.Time = (startDateTime:hours(1):endDateTime)'; % create time var
mr2.albedo(height(mr2.Time),1) = NaN; % create empty albedo column for data storage
merra2 = table2timetable(struct2table(mr2)); % convert to timetable
merra2.albedo(:,:) = NaN; % make sure that instead of 0 we have NaN


for i = 1:Nfiles % Loop for each nc-file

    fl = ncfiles(i).name; % extract the filename
    fprintf('Analyzing file %s\n', fl); % keep track of progress

    % Extract the year, month, and day from the filename
    fileDateStr = fl(28:35); % Extracts the date part from the filename
    fileDate = datetime(fileDateStr, 'InputFormat', 'yyyyMMdd'); % Convert to datetime
    
    % Check if the file date is within the range of the startDateTime and endDateTime
    if fileDate >= startDateTime && fileDate <= endDateTime
        % Calculate the index in merra2 for the file date
        index = find(merra2.Time == fileDate); % Find the index corresponding to the file date
        
        % Read albedo data from the file
        albedoData = ncread(fl, 'SNICEALB');

        % Place the albedo data into merra2 at the corresponding index
        merra2.albedo(index:3:index+23) = albedoData; % Adjust indices accordingly if needed
        
    else
        fprintf('File date %s is not within the specified range.\n', datestr(fileDate));
    end
end

clear albedoData fl i fileDate fileDateStr index mr2 startDateTime endDateTime

% % мой предыдущий код (снизу)
% 
% for i = 1:Nfiles % Loop for each nc-file
% 
%     fl       = ncfiles(i).name;     % extract the filename
%     ncInfo = ncinfo(fl); % get info on file
%     variableNames = {ncInfo.Variables.Name}; % get var names
%     
% 
%     fprintf('Analyzing file %s\n', fl); % keep track of progress
% 
% %     % Check if the file has all three variables
% %     if ~ismember('ALBEDO', variableNames) || ...
% %        ~ismember('TS', variableNames) || ...
% %        ~ismember('SWGDN', variableNames)
% %         fprintf('File %s is missing one or more variables.\n', fl);
% %         continue; % Skip to the next file
% %     end
% 
%     if ~ismember('SNICEALB',variableNames)
%         albedo = NaN;
%     else
%         albedo = ncread(fl,'SNICEALB');
%     end
% 
%     albedo   = ncread(fl,'SNICEALB'); % albedo
% 
%     % loop to extract hourly values from single .nc file
%     for nn = 1:files_time 
% %         indx = nn + (24*(i-1)); % which row
% %         mr2.albedo(indx,1) = albedo(1,1,nn); % put albedo values in mr2
% %         mr2.tsfc(indx,1)   = tsfc(1,1,nn);
% %         mr2.swd(indx,1)    = swd(1,1,nn);
%         % Calculate index for mr2
%         indx = (nn-1) * 3 + 1 + (24*(i-1)); % calculate index
%         
%         % Put albedo values in mr2
%         mr2.albedo(indx, 1) = albedo(1, 1, nn); % put albedo value
%         
%         % Fill NaN values for the 2 intermediate rows
%         mr2.albedo(indx + 1:indx + 2, 1) = NaN;
% 
%     end
% end



%% 3) separate years
% must extract only albedo and have a timetable for each year

% i) extract albedo
% mr2a = merra2(:,1);

% ii) Eradicate the mere notion of leap years 
feb29 = month(merra2.Time) == 2 & day(merra2.Time) == 29;
merra2 = merra2(~feb29,:); % update mr2a, now w/ no leap days

% iii) preparation for the year-separation loop 
yy = sy; % start year
endyy = ey+1; % to establish how many years total

% iv) find out the index of the row from which will start saving data
    % this is overly complex in case datasets don't start exactly on jan1,2009
date_to_find = datetime((yy), 1, 1); % specify start date
row_index = find(merra2.Properties.RowTimes == date_to_find); % bingo
% next year row
next_date = datetime((yy)+1, 1, 1); 
next_row = find(merra2.Properties.RowTimes == next_date);
next_row = 8761;
row_dist = next_row - row_index; % how many rows in full year


% c) separating the years * * * * * * * * * * * * * * * * * * * *
for i = 1:(endyy-sy) % 2009-2022 full years
    m = merra2(row_index:row_index+row_dist-1,:); % get all data for specified year
%     m = merra2;
    mx = m; % test variable

    % 1) find 1st existing value and fill all NaNs before
    fi = find(~isnan(m.albedo),1,'first'); % id the row of first existing value
    mx.albedo(1:fi) = m.albedo(fi); % fill 1st missing values with 1st valid value

    % 2) now fill the last missing values with the last valid value
    fi = find(~isnan(m.albedo),1,'last');
    mx.albedo(fi+1:height(m)) = m.albedo(fi);

    % 3) interpolate any missing data in between
    mx.albedo = fillmissing(mx.albedo,'linear');
%         m.albedo = fillmissing(m.albedo,'movmean',24); % alternative
    

    % 10) prep for next iteration
    row_index = row_index + row_dist; % update row_index
    nom = (['merra' num2str(yy)]); % specify file name
    yy = yy + 1; % update year
    assignin('base',char(nom),mx); % save separately
end

clear feb29 date_to_find row_index next_date next_row row_dist i m mx fi yy nom


%% Save the albedo so that IceModel v.1 can read it (added March 4th, 2024)
% Thus, rename each file to "Data", keep as timetable, save with naming
% convention of model_aws_year.mat (e.g., merra_kpcl_2013.mat).
% Model options: modis, merra, mar, racmo.
% Save to icemodel>input>userdata

% Before that, save figures of albedos just to visually inspect them.
% In figures>29>6

for yyyy=sy:ey
    fl = eval(['merra' num2str(yyyy)]); % get current file
    metVars = fl.Properties.VariableNames; % get var names
    tiempo = fl.Time; % the time
    FolderPath = [mp 'drive_python/figures/29. new aws/7. merra albs/']; % save figs here
    for v = 1:numel(metVars) % loop through each variable
        vnom = metVars{v}; % specify name of var about to select
        theVar = fl.(vnom); % extract data of the specified variable
        % begin plotting:
        figure('visible', 'off'); % create figure without displaying it; % create new fig for every var, year, and AWS
        plot(tiempo,theVar)
        xlabel('Time'); ylabel(vnom);
        title([aws ' ' vnom ' ' num2str(yyyy)])
        % save figure
        nom = [aws '_' vnom '_' num2str(yyyy) '.png']; % specify name
        FullPath = fullfile(FolderPath,nom); % establish full path
        print(FullPath,'-dpng',['-r' num2str(300)]); % save with 300 dpi
    end
end
clear fl metVars tiempo FolderPath v vnom theVar nom FullPath yyyy

% save each year's albedo to the folder that IceModel will enjoy
for i = sy:ey 
    nom = ['merra' num2str(i)]; % get name of file
    Data = eval(nom); % activate select file
    outloc = [mp 'input/userdata/']; % output location (save to this new place)
    newname = ['merra_' aws '_' num2str(i) '.mat'];
    nnal = [outloc newname]; % new name and location (for met file)
    save(nnal,'Data'); % save in the new location
end
clear i nom outloc newname nnal Data








%% old but maybe not entirely useless?
% %% just swgdn for just 2016
% 
% %% 1) get files ===========================================================
% % Read the names of nc-files in the folder 
% pth = [mp 'drive_python/input/MERRA2/M2T1NXRAD_5.12.4/kanl2016/*.nc'];
% ncfiles      = dir(pth);   % get basic info
% 
% Nfiles       = length(ncfiles);  % Obtain total number of .nc files 
% file_name    = ncfiles(1).name;  % extract 1st filename
% files_time   = length(ncread(file_name,'time')); % if 24 - hourly
% % ncdisp(ncfile.name); % to see ALL info of .nc file
% 
% 
% 
% %% 2) extract vars and put in mr2 structure ===============================
% % Nfiles = 3;
% for i = 1:Nfiles % Loop for each nc-file
% 
%     fl       = ncfiles(i).name;     % extract the filename
%     ncInfo = ncinfo(fl); % get info on file
%     variableNames = {ncInfo.Variables.Name}; % get var names
%     
% 
%     fprintf('Analyzing file %s\n', fl); % keep track of progress
% 
%     % check if albedo exists. if not, make NaN. if does, get values
%     if ~ismember('ALBEDO',variableNames)
%         albedo(1,1,24) = NaN;
%     else
%         albedo = ncread(fl,'ALBEDO');
%     end
% 
%     if ~ismember('TS',variableNames)
%         tsfc(1,1,24) = NaN;
%     else
%         tsfc = ncread(fl,'TS');
%     end
% 
%     if ~ismember('SWGDN',variableNames)
%         swd = NaN;
%     else
%         swd = ncread(fl,'SWGDN');
%     end
% 
% %     % Check if the file has all three variables
% %     if ~ismember('ALBEDO', variableNames) || ...
% %        ~ismember('TS', variableNames) || ...
% %        ~ismember('SWGDN', variableNames)
% %         fprintf('File %s is missing one or more variables.\n', fl);
% %         continue; % Skip to the next file
% %     end
% 
% 
% %     albedo   = ncread(fl,'ALBEDO'); % albedo
% %     tsfc     = ncread(fl,'TS');     % skin temp
% %     swd      = ncread(fl,'SWGDN');  % incoming shortwave radiation
% 
% 
%     % loop to extract hourly values from single .nc file
%     for nn = 1:files_time 
%         indx = nn + (24*(i-1)); % which row
%         mr2.albedo(indx,1) = albedo(1,1,nn); % put albedo values in mr2
%         mr2.tsfc(indx,1)   = tsfc(1,1,nn);
%         mr2.swd(indx,1)    = swd(1,1,nn);
% 
%     end
% end
% 
% % add time and convert
% % startDateTime = datetime(2009, 1, 1, 0, 0, 0);
% % endDateTime = datetime(2023, 1, 1, 23, 0, 0);
% % mr2.Time = (startDateTime:hours(1):endDateTime)';
% % merra2 = table2timetable(struct2table(mr2));