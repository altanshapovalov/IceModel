% Purpose: check new Greenland weather stations for quality (if gaps and plots)
% Created on Feb 21, 2024; Maxim Altan-Lu Shapovalov

% Contents:
% -1) read nc file
%  0) Get start and end time
%  1) extract and allocate data
%  2) further corrections
%  3) Separate Years and Save
%  4) Pit-stop – fix all NaNs and fill all gaps
%  5) Plot check for albedos
%  6) Compare parameters of kanm to kanl (may 29)


%% 
clear; clc
mp = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/'; % mp = max's path

% CHANGE THIS FOR DIFFERENT WEATHER STATION
% aws = 'KPC_L';1
% aws = 'KPC_U';2
% aws = 'NUK_L';3
% aws = 'NUK_U';4
% aws = 'QAS_L';5
% aws = 'QAS_M';6
% aws = 'QAS_U';7
% aws = 'SCO_L';8
% aws = 'SCO_U';9
% aws = 'THU_L';10
% aws = 'THU_U';11
% aws = 'UPE_L';12
% aws = 'UPE_U';13


awss = {'KPC_L' 'KPC_U' 'NUK_L' 'NUK_U' 'QAS_L' 'QAS_M' 'QAS_U' 'SCO_L' 'SCO_U'...
    'THU_L' 'THU_U' 'UPE_L' 'UPE_U'};
names = {'kpcl' 'kpcu' 'nukl' 'nuku' 'qasl' 'qasm' 'qasu' 'scol' 'scou' 'thul'...
    'thuu' 'upel' 'upeu'}; % for exporting met files

% awss = {'KANL' 'KANM' 'KANU'};
% names = {'kanl' 'kanm' 'kanu'};

% plotting colors:
color_codes = {...
    [0.298, 0.447, 0.690], ... % Dark blue
    [0.333, 0.659, 0.408], ... % Dark green
    [0.768, 0.305, 0.321], ... % Dark red
    [0.505, 0.447, 0.698], ... % Dark purple
    [0.800, 0.725, 0.454], ... % Dark yellow
    [0.392, 0.709, 0.803], ... % Dark cyan
    [0.690, 0.388, 0.474], ... % Dark magenta
    [0.220, 0.494, 0.717], ... % Dark cerulean
    [0.310, 0.749, 0.357], ... % Dark lime green
    [0.839, 0.372, 0.372], ... % Dark salmon
    [0.737, 0.741, 0.133], ... % Dark olive
    [0.282, 0.282, 0.282]  ... % Dark gray
    };


%% -1) read nc file ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

% CHANGE THIS: which KAN file
for xd = 9%:length(awss) % FOR NOW JUST TEST ON 1 AWS
    aws = awss{xd};
    sy = 2010; % start year
    ey = 2022; % end year
    
    if numel(awss)>3 % if not KAN, but the new ones
        pth = [mp 'drive_python/input/runoff_temp_feedback/others/' aws '_hour.nc'];
    else
        pth = [mp 'drive_python/input/runoff_temp_feedback/raw/' aws '_hour.nc'];
    end
    
    ncfile = dir(pth); % get basic info
    % ncdisp(ncfile.name); % to see ALL info of .nc file
    fnom = ncfile.name; % get file name
    % just in case previous line doesn't catch if file didn't read
    % if isempty(ncfile)
    %     warning(['Did not read ' aws ' file properly'])
    % else
    %     disp([aws ' file read successfully'])
    % end
    
    
    
    %% 0) Get start and end time –––––––––––––––––––––––––––––––––––––––––––––––––––
    
    fileInfo = ncinfo(fnom); % Get information about the file
    globalAttributes = fileInfo.Attributes; % Access the global attributes
    
    % Find the attribute values by name
    attributeNames = {'time_coverage_start', 'time_coverage_end'};
    attributeValues = cell(1, numel(attributeNames));
    
    % get the actual values
    for i = 1:numel(attributeNames)
        attributeName = attributeNames{i}; % have as one thing
        attributeValue = ''; 
    
        for j = 1:numel(globalAttributes)
            if strcmp(globalAttributes(j).Name, attributeName)
                attributeValue = globalAttributes(j).Value;
                break;
            end
        end
    
        attributeValues{i} = attributeValue;
    end
    
    % Display the attribute values
    for i = 1:numel(attributeNames)
        attributeName = attributeNames{i};
        attributeValue = attributeValues{i};
    
        if ~isempty(attributeValue)
    %         disp([attributeName, ': ', attributeValue]);
        else
            disp(['Attribute ', attributeName, ' not found.']);
        end
    end
    
    clear attributeName attributeNames attributeValue fileInfo globalAttributes i j
    
    
    
    %% 1) extract and allocate data ––––––––––––––––––––––––––––––––––––––––––––––––
    mett.tair = ncread(fnom,'t_u'); % air temp
    mett.swd = ncread(fnom,'dsr'); % shortwave down radiation
    mett.lwd = ncread(fnom,'dlr'); % longwave down
    mett.albedo = ncread(fnom,'albedo'); % albedo
    mett.snow(length(mett.albedo),:) = NaN;
    mett.rain(length(mett.albedo),:) = NaN;
    mett.melt(length(mett.albedo),:) = NaN;
    mett.runoff(length(mett.albedo),:) = NaN;
    mett.shf = ncread(fnom,'dshf_u'); % sensible heat flux
    mett.lhf = ncread(fnom,'dlhf_u'); % latent heat flux
    mett.smb(length(mett.albedo),:) = NaN;
    mett.snowd(length(mett.albedo),:) = NaN;
    mett.cfrac = ncread(fnom,'cc'); % cloud fraction
    mett.tsfc = ncread(fnom,'t_surf'); % surface temperature
    mett.psfc = ncread(fnom,'p_u'); % surface pressure
    mett.MODIS(length(mett.albedo),:) = NaN;
    mett.wspd = ncread(fnom,'wspd_u'); % wind speed
    mett.wdir = ncread(fnom,'wdir_u'); % wind direction
    mett.rh = ncread(fnom,'rh_u'); % relative humidity
    mett.date = ncread(fnom,'time'); % tiempo
    
    lat = ncread(fnom,'lat'); lon = ncread(fnom,'lon');
%     disp([aws ' ' num2str(lat) ' ' num2str(lon)]) % display coordinates
    
    
    
    
    %% 2) further corrections ––––––––––––––––––––––––––––––––––––––––––––––––––––––
    
    % get info on units from Matt's MAR
    load([mp 'input/matt_og_files/met/met_behar_MAR_2016_1hr.mat']);
    units   =   met.Properties.VariableUnits; clear met
    
    % Time management
    timeStartStr = attributeValues{1};  % time_coverage_start attribute value
    timeEndStr = attributeValues{2};    % time_coverage_end attribute value
    timeStart = datetime(timeStartStr); % , 'InputFormat', format);
    timeEnd = datetime(timeEndStr); %, 'InputFormat', format);
    % create time array
        T = (timeStart:hours(1):timeEnd)';
        % apply it to datafile (mett)
        mett.Time = T;
        met = table2timetable(struct2table(mett));
    
    % correct the NaN columns
    met.snow(met.snow==0) = NaN;
    met.rain(met.rain==0) = NaN;
    met.melt(met.melt==0) = NaN;
    met.runoff(met.runoff==0) = NaN;
    met.smb(met.smb==0) = NaN;
    met.snowd(met.snowd==0) = NaN;
    met.MODIS(met.MODIS==0) = NaN;
    
    % turn date into datenum (to match Matt's met)
    met.date = datenum(T);
    
    % convert some units
    met.tair = met.tair + 273.15; % Celcius to Kelvin
    met.tsfc = met.tsfc + 273.15;
    met.psfc = met.psfc .* 100; % hPa (hectopascals) -> Pa
    
    % update some properties
    met.Properties.VariableUnits = units; 
    
    % rename to free up the "met" name for the future
    full_met_08_23 = met; 
    
    clear T pth mett units ncfile fnom met attributeValues timeEnd timeEndStr ...
        timeStart timeStartStr lat lon
    
    
    
    
    %% 3) Separate Years and Save ––––––––––––––––––––––––––––––––––––––––––––––––––
    
    % Eradicate the mere notion of leap years
    feb29 = month(full_met_08_23.Time) == 2 & day(full_met_08_23.Time) == 29;
    full_met_noleap = full_met_08_23(~feb29,:);
    
    % здесь я сперва просто подготавливаю таблицу в которой буду содержать инфу о летних нанах
    % Создаём ячейчатую (cell) таблицу с именами колонок
    yrz = sy:2022; num_years = length(yrz); % количество лет
    endyy = 2023; % to establish how many years
    column_names = cell(1, num_years);  % Создаем пустую таблицу
    for year = sy:2022
        col_name = num2str(year);  % Преобразовываем год в строку
        column_names{year - (sy-1)} = col_name;  % Заполняем название колонки
    end
    % Добавляем имена рядов к этой таблице
    row_names = {'Total NaNs' 'Total NaNs (%)' 'Total Gaps' 'Max Gap Size' 'Avg Gap Size'};
%     row_names2 = {'Poor' 'Eh' 'Okay' 'Good' 'Perfect'};
    % Создаём пустую таблицу с указанными именами колонок и рядов
    data_table = cell(numel(row_names), num_years);
    naninfo_table = cell2table(data_table, 'VariableNames', column_names, 'RowNames', row_names);
%     data_table2 = cell(numel(row_names2), numel(awss));
%     awsinfo_table = cell2table(data_table2,'VariableNames',awss,'RowNames',row_names2);
    clear yrz num_years column_names col_name row_names data_table year % чистка

%         % __ preliminary stuff before the loop
%     % start year
%     yy = sy; % starting with 2009 (2010 for KANU) because first full year of data
%     endyy = 2023; % to establish how many years
%     % find out the index of the row from which will start saving data
%     date_to_find = datetime(yy, 1, 1); % specify start date
%     row_index = find(full_met_noleap.Properties.RowTimes == date_to_find); % bingo
%     % next year row
%     next_date = datetime(yy+1, 1, 1); 
%     next_row = find(full_met_noleap.Properties.RowTimes == next_date); 
%     row_dist = next_row - row_index; % how many rows in full year

    % Find the row index of the first instance of January 1st of any year
    jan_1st_indices = find(month(full_met_noleap.Time) == 1 & day(full_met_noleap.Time) == 1);
    % Get the row index of the first occurrence of January 1st
    row_index = jan_1st_indices(1);
    % Extract the year from the Time variable at the identified row index
    year_of_first_jan_1st = year(full_met_08_23.Time(row_index));
    yy = year_of_first_jan_1st; % have the bottom loop start w/ this year
    next_date = datetime(yy+1,1,1); % next year row
    next_row = find(full_met_noleap.Properties.RowTimes == next_date); 
    row_dist = next_row - row_index; % how many rows in full year
    % find what is the last full year
    year_of_last_value = year(full_met_08_23.Time(end)); % year of last value
    checkyear = year_of_last_value; % for testing purposes (can change)
    % see how many hours are in the last year of dataset
    checkkk = height(full_met_noleap(year(full_met_noleap.Time) == checkyear,:));
    if checkkk < 8760 % if less than full year
        ey = year_of_last_value - 1;
    end
    
    for i = 1:(ey+2 - yy) 
        m = full_met_noleap(row_index:row_index+row_dist-1,:); % get all data for specified year
    
        % IF YOU WANNA SEE NANS BEFORE PATCHING, UNCOMMENT THIS SECTION*********
%         % прежде чем что-либо залатывать, смотрим сколько НаНов
%         % создано 12 октября, сохраняю таблицу с инфой про летние наны)
%         % fyi, jja is june, july, august (so, summer)
%             % i) обозначь, с какми периодом времени мы работаем
% %         t1 = find(m.Time == datetime(yy,6,1,0,0,0)); % start of jja
% %         t2 = find(m.Time == datetime(yy,8,31,23,0,0)); % end of jja
%         t1 = find(m.Time == datetime(yy,4,1,0,0,0)); % start of jja
%         t2 = find(m.Time == datetime(yy,9,30,23,0,0)); % end of jja
%             % ii) вытащи только альбедо и вырежи только лето
%         albb = m.rh; sumalbb = albb(t1:t2,:);
%             % iii) расчёты
%         totalnans = sum(ismissing(sumalbb)); % integer of nans
%         porciento = (totalnans/height(sumalbb))*100; % процентное соотношение
%         % расчитываем прорехи
%             nan_lengths = zeros(size(sumalbb)); % приготовь документ в который пойдёт инфа о дырах
%             current_length = 0; % счётчик с нуля
%             for ii = 1:numel(sumalbb) % на всю длину
%                 if isnan(sumalbb(ii)) % если пункт нан..
%                     current_length = current_length + 1; % ..увеличь счётчик
%                 else % ежели не нан, а цифра, то..
%                     nan_lengths(ii) = current_length; % ..останавливай счётчик и сохряняй инфу
%                     current_length = 0; % счётчик сначала чтоб не смешивать дыры
%                 end
%             end
%         maxhole = max(nan_lengths); % самая большая/длинная дыра
%         totalholes = nan_lengths(nan_lengths > 0); % считай сколько дыр всецелом
%         avghole = mean(totalholes);
%         % conditional in case NaNs are for the entire jja
%         if totalnans == height(sumalbb)
%             totalholes = 1; % 1 big hole
%             maxhole = height(sumalbb); % entire height of sumalbb is the biggest hole
%             avghole = maxhole; % la misma
%         end
%         % another conditional in case no NaNs at all
%         if totalnans == 0
%             avghole = 0;
%         end
%             % iv) сохраняем в таблицу
%         naninfo_table{1,(yy-(sy-1))} = num2cell(totalnans);
%         naninfo_table{2,(yy-(sy-1))} = num2cell(round(porciento));
%         naninfo_table{3,(yy-(sy-1))} = num2cell(height(totalholes));
%         naninfo_table{4,(yy-(sy-1))} = num2cell(maxhole);
%         naninfo_table{5,(yy-(sy-1))} = num2cell(round(avghole));
        % END OF NAN INFO SECTION **********************************************
    
        % fix albedo by MGC's conventions *** *** *** *** ***
        % 1) find 1st existing value and fill all NaNs before
        fi = find(~isnan(m.albedo),1,'first'); % id the row of first existing value
        m.albedo(1:fi) = m.albedo(fi); % fill 1st missing values with 1st valid value
        % 2) now fill the last missing values with the last valid value
        fi = find(~isnan(m.albedo),1,'last');
        m.albedo(fi+1:height(m)) = m.albedo(fi);
        % 3) interpolate any missing data in between
%         m.albedo = fillmissing(m.albedo,'linear');
%         m.albedo = fillmissing(m.albedo,'movmean',24); % alternative
        % 4) set Nov, Dec, Jan, and Feb to 0.8 (bc Matt said so)
        febrow = find(m.Time == datetime(yy,3,1)) - 1; % last hr of feb
        novrow = find(m.Time == datetime(yy,11,1)); % first hr of nov
        endrow = height(m); % last hr of dec
        m.albedo(1:febrow) = 0.8; % make jan+feb 0.8
        m.albedo(novrow:endrow) = 0.8; % same with nov+dec
    
        % get rid of NaN rows for section 4 of the code: ++++ ++++ ++++
        m = removevars(m,{'snow','rain','melt','runoff','smb','snowd','MODIS','date'});
    
        row_index = row_index + row_dist;% update row_index
        nom = (['xmet' num2str(yy)]); % specify file name
        yy = yy + 1; % update year
        assignin('base',char(nom),m); % save separately

        checkname = ['xmet' num2str(ey)]; % name of last possible file
        if exist(checkname, 'var') % if last possible file exists
            break % break out of this loop and get to nan patching
        end

    end

%     disp(['NaN info for ' aws ':'])
%     disp(naninfo_table); 
    
    clear yy row_index row_dist nom next_row next_date m i feb29 date_to_find ...
        full_met_noleap full_met_08_23 fi endrow febrow novrow endyy maxhole...
        avghole totalholes t1 t2 nan_lengths current_length porciento totalnans...
        albb sumalbb ii row_names naninfo_table checkkk checkyear

    %% PATCHING NaN GAPS -------------------------------------------------------
    disp(['About to begin patching gaps for ' aws]) % track progress

    sy = year_of_first_jan_1st; % overwrite start year based on the earliest available full dataset
    for yy = sy:ey
        
        % get info
        fl = ['xmet' num2str(yy)]; % get the met name

        fl = eval(fl); % get the actual timetable info      
        nom = (['met' num2str(yy)]); % save as (in section 4d)
        % -----------------
        
        % check for NaNs and display it
        %     % ChatGPT method: Check for NaN values in each column
        %     nanCheck = varfun(@(x) any(ismissing(x)), fl); 
        %     disp(nanCheck); % Display the result
        
        % a) see how many NaNs there are exactly in each column ====================
        nanCount = sum(ismissing(fl)); % Count NaN values in each column
        % Create a table to display the result
        resultTable = table(nanCount', 'VariableNames', {'NaN_Count'});
        resultTable.Properties.RowNames = fl.Properties.VariableNames;
    %     disp(resultTable); % Display the result
        metVars = fl.Properties.VariableNames; % get var names
        sumnans = max(nanCount);% nans 
        window1 = 2160; % 3 months
        w2 = height(fl)-window1; w3 = height(fl);


        for v = 1:numel(metVars) % loop through each variable
            vnom = metVars{v}; % specify name of var about to select
            theVar = fl.(vnom); % extract data of the specified variable
            % begin patching
            % take care of first values for better patching methods
            if all(isnan(theVar(1:window1,:))) && all(isnan(theVar(w2:w3,:)))
                % if the 1 and last 3 months are empty, move on to next year
                disp([aws ' ' num2str(yy) ' ' vnom ' is bad, skipping year+var'])
            elseif ~isnan(theVar(1,:)) && ~isnan(theVar(end,:))
                % if both 1st and last values do exist, skip this loop entirely
            else % if not all is lost...
                % 1) find rows w/ the 1st and last existing data
                firstdatarow = find(~isnan(theVar),1,'first');
                lastdatarow = find(~isnan(theVar),1,'last');
                % 2) see which one is closest to its edge
                first_dist = firstdatarow - 1;
                last_dist = height(fl) - lastdatarow;
                if first_dist < last_dist || first_dist == last_dist
                    % if 1st value is closer to its edge than the last, OR if
                    % the two distances are equal, choose 1st value
                    winner = firstdatarow;
                    SearchArea = winner + 120; % look at the NEXT 5 days
                    subset = theVar(winner:SearchArea,:);
                    nonNaNsubset = subset(~isnan(subset)); % get rid of NaNs
                    bandaidNum = median(nonNaNsubset);
                elseif first_dist > last_dist
                    % if last value is closer to its edge, choose that
                    winner = lastdatarow;
                    SearchArea = winner - 120; % look at PREVIOUS 5 days
                    subset = theVar(SearchArea:winner,:);
                    nonNaNsubset = subset(~isnan(subset)); % get rid of NaNs
                    bandaidNum = median(nonNaNsubset);
                end
                % 3) fill missing 1st/last value (ensuring not to overwrite)
                if isnan(theVar(1,:)) % if 1st value is NaN
                    theVar(1,:) = bandaidNum; 
                end
                if isnan(theVar(end,:)) % if last value is NaN
                    theVar(height(fl),:) = bandaidNum;
                end
            end


            % Previous New * * * * * * * * * * * *
%             if isnan(theVar(1,:)) && any(~isnan(theVar(2:336,:))) % if the 1st
%                 % value is NaN but there is at least 1 non-nan value in the 1st
%                 % 336hrs (2weeks) of the dataset...
%                 row1 = find(~isnan(theVar(1:336)),1,'first'); % row with 1st non-nan data
%                 row2 = row1+120; % 120hrs (5days) in the future...
%                 bandaidNum = median(theVar(row1:row2,:)); % get median of existing numbers to use for 1st value
%                 theVar(1,:) = bandaidNum; % create first value based on the educated guess
%             elseif all(isnan(theVar(1:336,:))) % if 1st 2 weeks are NaN...
%                 lastrow1 = height(fl) - 120; % 5 days away from last day
%                 lastrow2 = height(fl); % the very last row
%                 bandaidNum = median(theVar(lastrow1:lastrow2,:)); % same as above
%                 theVar(1,:) = bandaidNum; % see purpose above
%             else
%                  % if first value exists, move on
%             end
%             % another conditional for when last doesn't exist but 1st does
%             if any(~isnan(theVar(1:720,:))) && all(isnan(theVar(8040:8760,:)))
%                 %if any numbers in 1st 30 days (720hrs) aren't NaN AND if all
%                 %the values in last months are NaN...
%                 row1 = find(~isnan(theVar(1:720)),1,'first'); % row with 1st non-nan data
%                 row2 = row1+120; % 120hrs (5days) in the future...
%                 lastrow2 = height(fl); % the very last row
%                 bandaidNum = median(theVar(row1:row2,:));
%                 theVar(lastrow2,:) = bandaidNum;
%             elseif all(isnan(theVar(1:720,:))) && all(isnan(theVar(8040:8760,:)))
%                 % if both 1st and last month are NaNs...
%                 disp([yy ' ' vnom ' for ' awss{xd} ' lacks 1st and last months. Will try to fix, but check if crud'])
%             end
%             % last check
%             if isnan(theVar(height(fl),:)) && any(~isnan(theVar(1:720,:))) % if 
%                 % the last value is NaN but there is at least 1 nonNan value in
%                 % the first month, take the median of that and make it the last
%                 row1 = find(~isnan(theVar(1:720)),1,'first'); % row with 1st non-nan data
%                 row2 = row1+120; % 120hrs (5days) in the future...
%                 bandaidNum = median(theVar(row1:row2,:)); % get median of existing numbers to use for 1st value
%                 theVar(height(fl),:) = bandaidNum; % create first value based on the educated guess
%             end

            % now that there's a solid 1st value, patch accordingly
            fl.(vnom) = fillmissing(theVar(:,:),'linear');
%             fl.(vnom) = fillmissing(theVar(:,:),'movmean',sumnans);

            
        end
        clear v vnom metVars theVar row1 bandaidNum lastrow1 lastrow2...
            firstdatarow lastdatarow first_dist last_dist winner SearchArea...
            subset nonNaNsubset
 
% used previously for patching gaps
%     
%         % b) fill missing values ===================================================
%         % max # of NaNs; this will be the window used for moving mean/median
%         mostnans = max(resultTable.NaN_Count);
%         
%         % if there are more than 90 NaNs, see how many NaNs are before/during summer
%         %       if none, then let it be
%         %       if some, patch em up using movmean or movmedian (compare)
%         sumrow = find(fl.Time == datetime(yy,9,1)); % id last summer day (sept 1)
%         
%         % set some cutoffs for deciding gaps and solution methods
%             biggap = 82; % one of the files had 81 nans, so i'm saying that anything 
%             % above 82 will be considered having not a tiny gap, but a big gap
%             tinygap = 12; % choosing 12 as cutoff because JR used window of 11 back in 
%             % one paper and i feel 6 hrs before + after is a reasonable window where 
%             % enough pattern can get covered
%         
%         if mostnans >= biggap
%             nanCount = sum(ismissing(fl((1:sumrow),:))); % Count NaN values in each column
%             % Create a table to display the result
%             resultTable = table(nanCount', 'VariableNames', {'NaN_Count pre-sept 1'});
%             resultTable.Properties.RowNames = fl.Properties.VariableNames;
%     %         disp(resultTable); % Display the result
%             sumnans = max(nanCount);% nans before+during summer
%         
%             if sumnans < tinygap % if barely any NaNs, linear
%                 fl(1:sumrow,:) = fillmissing(fl(1:sumrow,:),"linear"); 
%     %             disp('Pre-sept: TINY gap, linear interpolation applied')
%     %             disp(['Max NaN pre sept: ', num2str(max(sum(ismissing(fl((1:sumrow),:)))))]);
%     %             disp(['Max NaN post sommeil: ', num2str(max(sum(ismissing(fl((sumrow:end),:)))))]);
%         
%         
%             elseif sumnans >= tinygap && sumnans <= biggap % if medium #, movmean
%                 fl(1:sumrow,:) = fillmissing(fl((1:sumrow),:),'movmean',hours(sumnans));
%     %             disp('Pre-sept: MEDIUM gap, mean interpolation applied')
%     %             disp(['Max NaN pre sept: ', num2str(max(sum(ismissing(fl((1:sumrow),:)))))]);
%     %             disp(['Max NaN post sommeil: ', num2str(max(sum(ismissing(fl((sumrow:end),:)))))]);
%         
%             elseif sumnans > biggap % if big gap, still use mean
%         %         fl(1:sumrow,:) = fillmissing(fl,'movmedian',hours(sumnans)); % nah...
%                 fl(1:sumrow,:) = fillmissing(fl((1:sumrow),:),'movmean',hours(sumnans));
% %                 fl(1:sumrow,:) = fillmissing(fl((1:sumrow),:),'linear');
%     %             disp('Pre-sept: BIG gap, mean interpolation applied')
%     %             disp(['Max NaN pre sept: ', num2str(max(sum(ismissing(fl((1:sumrow),:)))))]);
%     %             disp(['Max NaN post sommeil: ', num2str(max(sum(ismissing(fl((sumrow:end),:)))))]);
%         
%         
%                 % alt:
%         %         fl = fillmissing(fl,'movmean',hours(sumnans));
%         %         disp('Pre-sept: BIG gap, mean interpolation applied')
%             elseif sumnans == 0 % if no NaNs in summer
%                 % move on to lower conditional
%             end
%         
%         % For smaller gaps:
%         elseif mostnans >= tinygap && mostnans < biggap
%             fl = fillmissing(fl,'movmean',hours(mostnans));    
%             % fl = fillmissing(fl,'movmedian',hours(mostnans)); % alternative
%     %         disp('ALL: MEDIUM gap, mean was used') % confirm
%         
%         elseif mostnans < tinygap
%             fl = fillmissing(fl,'linear'); 
%     %         disp('ALL: TINY gap, linear was used') 
%         
%         end
%         
%         % fix gap after summer
%         PostSummerGap = num2str(max(sum(ismissing(fl((sumrow:end),:)))));
%         if PostSummerGap > 0
%             gapnan = max(sum(ismissing(fl))); % Count NaN values in each column
%             fl = fillmissing(fl,"movmean",hours(gapnan+1));
%             % check for gaps
%             gapnan = max(sum(ismissing(fl))); % check again
%             count = 0; % keep count just for fun of the visualization
%         
%             while gapnan > 0 % run this until gapnan becomes zero
%                 fl = fillmissing(fl,"movmean",hours(gapnan+1)); % patch again
%                 gapnan = max(sum(ismissing(fl))); % check again
%                 count = count + 1; % count iterations
%     %             disp(['Iteration ' num2str(count) ', Current value of gapnan: ', num2str(gapnan)]);
%             end
%         
%             if gapnan == 0 % final check if there are any NaNs left
%     %             disp('Post-summer patched; please check')
%             end
%         
%         else
%     %         disp('No NaN gaps post-summer; please double-check')
%         end
%     
        nantest = num2str(max(sum(ismissing(fl(:,:)))));
        disp([nantest ' NaNs in ' num2str(yy)])
        assignin('base',char(nom),fl); % save separately
        
    end

    % clear the xmet timetables (чтобы не мешались)
    for year = 2008:2023
        variable_name = ['xmet', num2str(year)];
        if exist(variable_name, 'var')
            clear(variable_name);
        end
    end

    clear nanCount resultTable mostnans ans biggap count gapnan PostSummerGap...
        sumnans sumrow tinygap yy nom nantest fl maxhole avghole totalholes t1 t2...
        nan_lengths current_length porciento totalnans albb sumalbb i year...
        variable_name


    %% Plot each var and save (for visual inspection)
    disp(['Creating plots for ' aws])

    for yyyy=sy:ey
        fl = eval(['met' num2str(yyyy)]); % get current file
        metVars = fl.Properties.VariableNames; % get var names
        tiempo = fl.Time; % the time
        FolderPath = [mp 'drive_python/figures/29. new aws/2. var check/']; % save figs here
        for v = 1:numel(metVars) % loop through each variable
            vnom = metVars{v}; % specify name of var about to select
            theVar = fl.(vnom); % extract data of the specified variable
            % begin plotting:
            figure('visible', 'off'); % create figure without displaying it; % create new fig for every var, year, and AWS
            plot(tiempo,theVar,'Color',color_codes{v})
            xlabel('Time'); ylabel(vnom);
            title([aws ' ' vnom ' ' num2str(yyyy)])
            % save figure
            nom = [aws '_' vnom '_' num2str(yyyy) '.png']; % specify name
            FullPath = fullfile(FolderPath,nom); % establish full path
            print(FullPath,'-dpng',['-r' num2str(300)]); % save with 300 dpi
        end
        close all %close figures once all are printed and saved
    end
    close all %close figures once all are printed and saved

    clear fl metVars tiempo FolderPath v vnom theVar nom FullPath

    % save each met file
    for god = sy:ey
        met = eval(['met' num2str(god)]);
        outloc = [mp 'input/met/'];
        newname = ['met_' names{xd} '_' names{xd} '_' num2str(god) '_1hr.mat'];
        nnal = [outloc newname]; % new name and location (for met file)
        save(nnal,'met'); % save in the new location
    end

    clear god met outloc newname nnal year_of_first_jan_1st

    disp(['Finished work on ' aws]) % keep track of progress

end

