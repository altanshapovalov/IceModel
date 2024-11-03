% New MODIS Albedo – extracting and organizing the data
% Created on May 16, 2024, by Maxim Altan-Lu Shapovalov

% Purpose: take the modis albedo extracted via Google Earth Engine Code Editor
% and organize 5 and 10 km albedo data for each AWS. Fyi, using the MCD43A3 
% MODIS product – Black Sky Shortwave (bc "broadband" and used by Alexander 
% et al., 2014)

% Ran on May 17th, 2024. Saved in icemodel>input>userdata>mcd43a3 (albedo data)

% Содержание:
% 0) Преамбула
% 1) Import csv file
% 2) Extracting albedo and establishing time
% 3) Expand & separate
% 4) Save the albedo so that IceModel v.1 can read it 


%% Преамбула
clear; clc
mp = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/'; % mp = max's path

FolderPath = [mp 'drive_python/figures/29. new aws/34. mcd43a3 single pixel/1. albedo/']; % save figs here

% CHANGE THIS FOR DIFFERENT WEATHER STATION
aws = 'kanl';
% aws = 'kanm';
% aws = 'kanu';
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
mets = {'kanl','kanm','kanu','kpcl','kpcu','nukl','nuku','qasl','qasu','scol',...
    'scou', 'thul', 'thuu', 'upel', 'upeu'};





%% 1) Import csv file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for xx = 1:length(mets)
    aws = mets{xx};
    % import
    datatable = readtable([mp 'drive_python/input/MODIS/MCD43A3/' aws '.csv'],...
        'ReadRowNames',true); 
    
    % establish start and end years based on what's inside the dataset
    % for start year (sy), find 1st instance of -01-01 in a cell
    rowNames = datatable.Properties.RowNames; % Extract row names from the table
    % Find the first instance of jan 1st
    for i = 1:length(rowNames)
        if contains(rowNames{i}, '_01_01')
            sy = str2double(rowNames{i}(1:4)); % Extract the first 4 numbers
            break; % Exit the loop once found
        end
    end
    % Find the last instance of dec 31st and have that as end year (ey)
    for i = length(rowNames):-1:1
        if contains(rowNames{i}, '_12_31')
            ey = (str2double(rowNames{i}(1:4))); % Extract the first 4 numbers
            break; % Exit the loop once found
        end
    end
    clear i rowNames
    
    
    %% 2) Extracting albedo and establishing time ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    T = (datetime(sy,1,1):days(1):datetime(ey,12,31))'; % create fitting time
    
    % PROBLEM: It appears that the modis dataset has some missing days, so the
        % rest of this section (2) is making sure that ALL days between 2009
        % and 2022 are present and those days have NaN for albedo
    
    % figuring out missing dates vvvvvvvv
    
    % ensure unique dates in each
    uq_dt = unique(datatable(:,1)); % dates from datatable
    uq_mod = unique(T); % dates from created time array (with full dates)
    dtt = table2array(uq_dt); % convert
    dt = datetime(dtt(:,1), 'InputFormat', 'yyyy_MM_dd'); % convert again
    
    alb10km = table2array(datatable(:,6)); % extract 10km albedo and convert to cell
    % alb5km  = table2array(datatable(:,2)); % extract 5km albedo and convert to cell
    
    % have a table with modis dates and albedo; still few dates missing
    data10 = table(dt,alb10km,'VariableNames',{'Time', 'albedo'});
    % data5 = table(dt,alb5km,'VariableNames',{'Time', 'albedo'});
    
    missingDates = setdiff(uq_mod, dt); % find out which dates are missing
    % Create table with missing dates and NaN albedo values
    missingData = table(missingDates, NaN(size(missingDates)), 'VariableNames', {'Time', 'albedo'});
    
    % Concatenate the missing data with the original dataset and put in
        % chronological order
    data10_upd = sortrows([data10; missingData], 'Time');
    % data5_upd = sortrows([data5; missingData], 'Time');
    
    mod10 = table2timetable(data10_upd); % convert to timetable
    % mod5 = table2timetable(data5_upd); % convert to timetable
    mod10.albedo = mod10.albedo / 1000; % put into albedo units
    % mod5.albedo = mod5.albedo / 1000; % put into albedo units
    
    % Get rid of data that is beyond my end year (ey)
    years10 = year(mod10.Properties.RowTimes); % Extract years from timetable dates
    % years5 = year(mod5.Properties.RowTimes);
    rowsToKeep10 = years10 <= ey; % Find rows where the year is above the threshold
    % rowsToKeep5  = years5 <= ey; % Find rows where the year is above the threshold
    filteredTable10 = mod10(rowsToKeep10,:); % Filter the timetable based on the rows to keep
    % filteredTable5 = mod5(rowsToKeep5,:); % Filter the timetable based on the rows to keep
    mod10 = filteredTable10; % update mod
    % mod5 = filteredTable5; % update mod
    
    % SUCCESS: now "mod" is a timetable that has all days and all albedo
        % values. Next: separate years, expand from daily to hourly, fill in
        % NaN gaps, and apply necessary corrections
    
    clear alb10km alb5km data10 data5 datatable dt dtt missingData missingDates T uq_dt...
        uq_mod data10_upd data5_upd years10 years5 rowsToKeep10 rowsToKeep5...
        filteredTable10 filteredTable5
    
    
    %% 3) Expand & separate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    for ixi = 1%:2 % efficient work of both 10 and 5 km albedo data
        if ixi == 1
            mod = mod10; % 1st work with mod10
        elseif ixi == 2
            mod = mod5; % then with 5km
        end
        % a) Eradicate the mere notion of leap years  * * * * * * * * * * 
        feb29 = month(mod10.Time) == 2 & day(mod10.Time) == 29;
        mod10 = mod10(~feb29,:); % update mod, now w/ no leap days
    %     mod5 = mod5(~feb29,:); % update mod, now w/ no leap days
        
        % b) preparation for the year-separation loop * * * * * * * * * *
        yy = sy; % start year
        endyy = ey+1; % to establish how many years total
        
        % find out the index of the row from which will start saving data
            % this is overly complex in case datasets don't start exactly on jan1,2009
        date_to_find = datetime(yy, 1, 1); % specify start date
        row_index = find(mod10.Properties.RowTimes == date_to_find); % bingo
        
        % next year row
        next_date = datetime(yy+1, 1, 1); 
        next_row = find(mod10.Properties.RowTimes == next_date); 
        row_dist = next_row - row_index; % how many rows in full year
        
        % set upper and lower limits of albedo for the loop
        uplim = 0.84;
        lowlim = 0.30;
        % NOTE: upper threshold chosen based on Alexander et al. 2014 (ctrl+F 0.84)
        
        % c) separating the years * * * * * * * * * * * * * * * * * * * *
    
        % separate the years/manipulate as necessary:
        for i = 1:(endyy-yy) 
            m = mod(row_index:row_index+row_dist-1,:); % get all data for specified year
    %         mx = m; % test variable
        
            % 1) find 1st existing value and fill all NaNs before + + + 
            fi = find(~isnan(m.albedo),1,'first'); % id the row of first existing value
            % check if very 1st value is NaN; if yes, make upper limit (0.84)
            if isnan(m.albedo(1,1))
                m.albedo(1,1) = uplim;
            end
            % fill gap via linear interp for a smoother transition
            m.albedo(1:fi) = fillmissing(m.albedo(1:fi),'linear');
        
            % 2) now fill the last missing values with the last valid value + + + 
            fi = find(~isnan(m.albedo),1,'last');
            % check if very last value is NaN; if yes, make upper limit (0.84)
            if isnan(m.albedo(end,1))
                m.albedo(end,1) = uplim;
            end
            % fill gap via linear interp for a smoother transition
            m.albedo(fi:end) = fillmissing(m.albedo(fi:end),'linear');
    
            % 3) interpolate any missing data in between + + +
            m.albedo = fillmissing(m.albedo,'linear');
    
            % 4) apply upper/lower limits
            m.albedo(m.albedo > uplim) = uplim; % upper
            m.albedo(m.albedo < lowlim) = lowlim; % lower
            m.albedo = movmean(m.albedo,11); % add movmean to daily dataset to exclude gnarly outliers (if any)
    
            % 5) expand (daily -> hourly)
            m1h_x = retime(m,'hourly','linear');
            T = m1h_x.Time;
            if length(T) > 8760 % if leap year
                feb29 = month(T) == 2 & day(T) == 29;
                m1h_x = m1h_x(~feb29,:); % find feb 29 and expunge it
            end
            % m1h is missing 23 rows/hours so must add time/data
            % 5.1) extract albedo and fill it up
            albedo = m1h_x.albedo;
            albedo(end:end+22,1) = NaN; % add next 22 hrs as NaN
            albedo(end+1,1) = uplim; % make last row uplim (0.84);
            albedo = fillmissing(albedo,'Linear'); % patch NaN gap
            % CHANGE A FEW MONTHS TO 0.84:
            % why: mcd43a3 has some odd troughs in earlier months when it
            % doesn't make sense for ice to be exposed. Thus, will constrain
            % albedo at colder months
            apr1row = 2161;
            oct1row = 6553; % just searched up rows manually for time's sake
            albedo(1:apr1row,1) = 0.84;
            albedo(oct1row:end,1) = 0.84;
            % 5.2) create a separate time var and make it leap-year-free
            myear = year(m1h_x.Time(1));
            Time = (datetime(myear,1,1,0,0,0):hours(1):datetime(myear,12,31,23,0,0))';
            if length(Time) > 8760 % if leap year
                feb29 = month(Time) == 2 & day(Time) == 29;
                Time = Time(~feb29,:); % find feb 29 and expunge it
            end
            % 5.3) convert to timetable
            m1h = timetable(Time,albedo);
        
            % 6) prep for next iteration
            row_index = row_index + row_dist; % update row_index
            if ixi == 1
                nom = (['modis_' num2str(yy)]); % specify file name
            else
    %             nom = (['modis_5km_' num2str(yy)]); % specify file name
            end
            yy = yy + 1; % update year
            assignin('base',char(nom),m1h); % save separately
        end
    
    end
    
    
    clear feb29 yy date_to_find row_index next_row next_date row_dist i m ...
        row_index nom endrow febrow fi hrs indx mdss mds mod modis_alb...
        modis_row mx Nfiles novrow T thresh thresh_min modis minn maxx endyy...
        ixi myear Time uplim lowlim albedo m1h m1h_x mod5 mod10
    
    
    %% 4) Save the albedo so that IceModel v.1 can read it ~~~~~~~~~~~~~~~~~~~~~~~~~
    % Thus, rename each file to "Data", keep as timetable, save with naming
    % convention of model_aws_year.mat (e.g., modis_kpcl_2013.mat).
    % Model options: modis, merra, mar, racmo.
    % Save to icemodel>input>userdata
    
    % Before that, save figures of albedos just to visually inspect them.
    % In figures>29>27
    disp('Saving figure')
    for yyyy=sy:ey
        fl = eval(['modis_' num2str(yyyy)]); % get current file
        metVars = fl.Properties.VariableNames; % get var names
        tiempo = fl.Time; % the time
        for v = 1:numel(metVars) % loop through each variable
            vnom = metVars{v}; % specify name of var about to select
            theVar = fl.(vnom); % extract data of the specified variable
            % begin plotting:
            figure('visible', 'off'); % create figure without displaying it; % create new fig for every var, year, and AWS
            plot(tiempo,theVar)
            xlabel('Time'); ylabel(vnom);
            title([aws ' ' vnom ' ' num2str(yyyy)])
            % save figure
            nom = [aws '_' num2str(yyyy) '.png']; % specify name
            FullPath = fullfile(FolderPath,nom); % establish full path
            print(FullPath,'-dpng',['-r' num2str(300)]); % save with 300 dpi
        end
    end
    clear fl metVars tiempo v vnom theVar nom FullPath yyyy
    
    % save each year's albedo to the folder that IceModel will enjoy
    for i = sy:ey 
        nom = ['modis_' num2str(i)]; % get name of file
        Data = eval(nom); % activate select file
        outloc = [mp 'input/userdata/mcd43a3/single pixel/']; % output location (save to this new place)
        newname = ['modis_' aws '_' num2str(i) '.mat'];
        nnal = [outloc newname]; % new name and location (for met file)
        save(nnal,'Data'); % save in the new location
    end
    clear i nom outloc newname nnal Data
    
    disp([aws ' completed'])
    disp(' ')

end





