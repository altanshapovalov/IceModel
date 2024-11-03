%% Read PROMICE AWS Data

% Purpose – get data from .nc files of KAN M/L/U automatic weather
%   stations. Specifically, must get hourly data of air temperature,
%   albedo, and overall prepare the data so that it can be forced into
%   Matt's IceModel.

% Started on April 26th, 2023

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
aws = 'KANL'; % options: KANM, KANL, KANU

sy = 2010; % start year

% % establish the start year (sy)
% if strcmp(aws,'KANM') || strcmp(aws,'KANL') % KANL and KANM start with 2009
%     sy = 2010; % было 2009, поменял на 2010 ибо в манускрипте всё равно на 2009 не смотрю
% elseif strcmp(aws,'KANU') % KANU starts with 2010
%     sy = 2010;
% else
%     error([aws ' is NOT a valid weather station on the K-transect (options: KANL, KANM, KANU)']);
% end




%% -1) read nc file ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

% CHANGE THIS: which KAN file
pth = [mp 'drive_python/input/runoff_temp_feedback/raw/' aws '_hour.nc'];

ncfile = dir(pth); % get basic info
ncdisp(ncfile.name); % to see ALL info of .nc file
fnom = ncfile.name; % get file name





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
        disp([attributeName, ': ', attributeValue]);
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





%% 2) further corrections ––––––––––––––––––––––––––––––––––––––––––––––––––––––

% get info on units from Matt's MAR
load([mp 'input/met/met_behar_MAR_2016_1hr.mat']);
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

%% ***** ***** SAVED AS:

load([mp 'drive_python/outputs/runoff_temp_feedback/' aws '/read_Promice_AWS/full_met_08_23.mat']); 





%% 3) Separate Years and Save ––––––––––––––––––––––––––––––––––––––––––––––––––

% Eradicate the mere notion of leap years
feb29 = month(full_met_08_23.Time) == 2 & day(full_met_08_23.Time) == 29;
full_met_noleap = full_met_08_23(~feb29,:);

% __ preliminary stuff before the loop
% start year
yy = sy; % starting with 2009 (2010 for KANU) because first full year of data
endyy = 2023; % to establish how many years
% find out the index of the row from which will start saving data
date_to_find = datetime(yy, 1, 1); % specify start date
row_index = find(full_met_noleap.Properties.RowTimes == date_to_find); % bingo
% next year row
next_date = datetime(yy+1, 1, 1); 
next_row = find(full_met_noleap.Properties.RowTimes == next_date); 
row_dist = next_row - row_index; % how many rows in full year


% здесь я сперва просто подготавливаю таблицу в которой буду содержать инфу о летних нанах
% Создаём ячейчатую (cell) таблицу с именами колонок
yrz = sy:2022; num_years = length(yrz); % количество лет
column_names = cell(1, num_years);  % Создаем пустую таблицу
for year = sy:2022
    col_name = num2str(year);  % Преобразовываем год в строку
    column_names{year - (sy-1)} = col_name;  % Заполняем название колонки
end
% Добавляем имена рядов к этой таблице
row_names = {'Total NaNs' 'Total NaNs (%)' 'Total Gaps' 'Max Gap Size' 'Avg Gap Size'};
% Создаём пустую таблицу с указанными именами колонок и рядов
data_table = cell(numel(row_names), num_years);
naninfo_table = cell2table(data_table, 'VariableNames', column_names, 'RowNames', row_names);
clear yrz num_years column_names col_name row_names data_table % чистка


for i = 1:(endyy - yy) % 2009-2022 full years
    m = full_met_noleap(row_index:row_index+row_dist-1,:); % get all data for specified year

    % прежде чем что-либо залатывать, смотрим сколько НаНов
    % создано 12 октября, сохраняю таблицу с инфой про летние наны)
    % fyi, jja is june, july, august (so, summer)
        % i) обозначь, с какми периодом времени мы работаем
    t1 = find(m.Time == datetime(yy,6,1,0,0,0)); % start of jja
    t2 = find(m.Time == datetime(yy,8,31,23,0,0)); % end of jja
        % ii) вытащи только альбедо и вырежи только лето
    albb = m.albedo; sumalbb = albb(t1:t2,:);
        % iii) расчёты
    totalnans = sum(ismissing(sumalbb)); % integer of nans
    porciento = (totalnans/height(sumalbb))*100; % процентное соотношение
    % расчитываем прорехи
        nan_lengths = zeros(size(sumalbb)); % приготовь документ в который пойдёт инфа о дырах
        current_length = 0; % счётчик с нуля
        for ii = 1:numel(sumalbb) % на всю длину
            if isnan(sumalbb(ii)) % если пункт нан..
                current_length = current_length + 1; % ..увеличь счётчик
            else % ежели не нан, а цифра, то..
                nan_lengths(ii) = current_length; % ..останавливай счётчик и сохряняй инфу
                current_length = 0; % счётчик сначала чтоб не смешивать дыры
            end
        end
    maxhole = max(nan_lengths); % самая большая/длинная дыра
    totalholes = nan_lengths(nan_lengths > 0); % считай сколько дыр всецелом
    avghole = mean(totalholes);
    % conditional in case NaNs are for the entire jja
    if totalnans == height(sumalbb)
        totalholes = 1; % 1 big hole
        maxhole = height(sumalbb); % entire height of sumalbb is the biggest hole
        avghole = maxhole; % la misma
    end
    % another conditional in case no NaNs at all
    if totalnans == 0
        avghole = 0;
    end
        % iv) сохраняем в таблицу
    naninfo_table{1,(yy-(sy-1))} = num2cell(totalnans);
    naninfo_table{2,(yy-(sy-1))} = num2cell(round(porciento));
    naninfo_table{3,(yy-(sy-1))} = num2cell(height(totalholes));
    naninfo_table{4,(yy-(sy-1))} = num2cell(maxhole);
    naninfo_table{5,(yy-(sy-1))} = num2cell(round(avghole));


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

        % done *** *** *** *** ***

            % get rid of NaN rows for section 4 of the code: ++++ ++++ ++++
            % snow, rain, melt, runoff, smb, snowd, MODIS
            m = removevars(m,{'snow','rain','melt','runoff','smb','snowd','MODIS'});
            % Save as "xmet" – done ++++ ++++ ++++


    row_index = row_index + row_dist;% update row_index
    nom = (['xmet' num2str(yy)]); % specify file name
    yy = yy + 1; % update year
    assignin('base',char(nom),m); % save separately
end

clear yy row_index row_dist nom next_row next_date m i feb29 date_to_find ...
    full_met_noleap full_met_08_23 fi endrow febrow novrow endyy maxhole...
    avghole totalholes t1 t2 nan_lengths current_length porciento totalnans...
    albb sumalbb ii




%% 4) Pit-stop – fix all NaNs and fill all gaps –––––––––––––––––––––––––––––––– (started may 15, 2023)

% load in the xmets with which I shall play and discover where NaNs are/n't

load([mp 'drive_python/outputs/runoff_temp_feedback/' aws '/read_Promice_AWS/xmets_2009_2022.mat']); %


%% 4 cont.

% % здесь я сперва просто подготавливаю таблицу в которой буду содержать инфу о летних нанах
% % Создаём ячейчатую (cell) таблицу с именами колонок
% yrz = sy:2022; num_years = length(yrz); % количество лет
% column_names = cell(1, num_years);  % Создаем пустую таблицу
% for year = sy:2022
%     col_name = num2str(year);  % Преобразовываем год в строку
%     column_names{year - (sy-1)} = col_name;  % Заполняем название колонки
% end
% % Добавляем имена рядов к этой таблице
% row_names = {'Total NaNs' 'Total NaNs (%)' 'Total Gaps' 'Max Gap Size' 'Avg Gap Size'};
% % Создаём пустую таблицу с указанными именами колонок и рядов
% data_table = cell(numel(row_names), num_years);
% naninfo_table = cell2table(data_table, 'VariableNames', column_names, 'RowNames', row_names);
% clear yrz num_years column_names col_name row_names data_table % чистка


for yy = sy:2022
    
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

%     % a_extra (создано 12 октября, сохраняю таблицу с инфой про летние наны)====
%     % fyi, jja is june, july, august (so, summer)
%         % i) обозначь, с какми периодом времени мы работаем
%     t1 = find(fl.Time == datetime(yy,6,1,0,0,0)); % start of jja
%     t2 = find(fl.Time == datetime(yy,8,31,23,0,0)); % end of jja
%         % ii) вытащи только альбедо и вырежи только лето
%     albb = fl.albedo; sumalbb = albb(t1:t2,:);
%         % iii) расчёты
%     totalnans = sum(ismissing(sumalbb)); % integer of nans
%     porciento = (totalnans/height(sumalbb))*100; % процентное соотношение
%     % расчитываем прорехи
%         nan_lengths = zeros(size(sumalbb)); % приготовь документ в который пойдёт инфа о дырах
%         current_length = 0; % счётчик с нуля
%         for i = 1:numel(sumalbb) % на всю длину
%             if isnan(sumalbb(i)) % если пункт нан..
%                 current_length = current_length + 1; % ..увеличь счётчик
%             else % ежели не нан, а цифра, то..
%                 nan_lengths(i) = current_length; % ..останавливай счётчик и сохряняй инфу
%                 current_length = 0; % счётчик сначала чтоб не смешивать дыры
%             end
%         end
%     maxhole = max(nan_lengths); % самая большая/длинная дыра
%     totalholes = nan_lengths(nan_lengths > 0); % считай сколько дыр всецелом
%     avghole = mean(totalholes);
%     % conditional in case NaNs are for the entire jja
%     if totalnans == height(sumalbb)
%         totalholes = 1; % 1 big hole
%         maxhole = height(sumalbb); % entire height of sumalbb is the biggest hole
%         avghole = maxhole; % la misma
%     end
%     % another conditional in case no NaNs at all
%     if totalnans == 0
%         avghole = 0;
%     end
%     % сохраняем в таблицу
%     naninfo_table{1,(yy-(sy-1))} = num2cell(totalnans);
%     naninfo_table{2,(yy-(sy-1))} = num2cell(porciento);
%     naninfo_table{3,(yy-(sy-1))} = num2cell(height(totalholes));
%     naninfo_table{4,(yy-(sy-1))} = num2cell(maxhole);
%     naninfo_table{5,(yy-(sy-1))} = num2cell(avghole);
    
    
    % b) fill missing values ===================================================
    % max # of NaNs; this will be the window used for moving mean/median
    mostnans = max(resultTable.NaN_Count);
    
    % if there are more than 90 NaNs, see how many NaNs are before/during summer
    %       if none, then let it be
    %       if some, patch em up using movmean or movmedian (compare)
    sumrow = find(fl.Time == datetime(yy,9,1)); % id last summer day (sept 1)
    
    % set some cutoffs for deciding gaps and solution methods
        biggap = 82; % one of the files had 81 nans, so i'm saying that anything 
        % above 82 will be considered having not a tiny gap, but a big gap
        tinygap = 12; % choosing 12 as cutoff because JR used window of 11 back in 
        % one paper and i feel 6 hrs before + after is a reasonable window where 
        % enough pattern can get covered
    
    if mostnans >= biggap
        nanCount = sum(ismissing(fl((1:sumrow),:))); % Count NaN values in each column
        % Create a table to display the result
        resultTable = table(nanCount', 'VariableNames', {'NaN_Count pre-sept 1'});
        resultTable.Properties.RowNames = fl.Properties.VariableNames;
%         disp(resultTable); % Display the result
        sumnans = max(nanCount);% nans before+during summer
    
        if sumnans < tinygap % if barely any NaNs, linear
            fl(1:sumrow,:) = fillmissing(fl(1:sumrow,:),"linear"); 
%             disp('Pre-sept: TINY gap, linear interpolation applied')
%             disp(['Max NaN pre sept: ', num2str(max(sum(ismissing(fl((1:sumrow),:)))))]);
%             disp(['Max NaN post sommeil: ', num2str(max(sum(ismissing(fl((sumrow:end),:)))))]);
    
    
        elseif sumnans >= tinygap && sumnans <= biggap % if medium #, movmean
            fl(1:sumrow,:) = fillmissing(fl((1:sumrow),:),'movmean',hours(sumnans));
%             disp('Pre-sept: MEDIUM gap, mean interpolation applied')
%             disp(['Max NaN pre sept: ', num2str(max(sum(ismissing(fl((1:sumrow),:)))))]);
%             disp(['Max NaN post sommeil: ', num2str(max(sum(ismissing(fl((sumrow:end),:)))))]);
    
        elseif sumnans > biggap % if big gap, still use mean
    %         fl(1:sumrow,:) = fillmissing(fl,'movmedian',hours(sumnans)); % nah...
            fl(1:sumrow,:) = fillmissing(fl((1:sumrow),:),'movmean',hours(sumnans));
%             disp('Pre-sept: BIG gap, mean interpolation applied')
%             disp(['Max NaN pre sept: ', num2str(max(sum(ismissing(fl((1:sumrow),:)))))]);
%             disp(['Max NaN post sommeil: ', num2str(max(sum(ismissing(fl((sumrow:end),:)))))]);
    
    
            % alt:
    %         fl = fillmissing(fl,'movmean',hours(sumnans));
    %         disp('Pre-sept: BIG gap, mean interpolation applied')
        elseif sumnans == 0 % if no NaNs in summer
            % move on to lower conditional
        end
    
    % For smaller gaps:
    elseif mostnans >= tinygap && mostnans < biggap
        fl = fillmissing(fl,'movmean',hours(mostnans));    
        % fl = fillmissing(fl,'movmedian',hours(mostnans)); % alternative
%         disp('ALL: MEDIUM gap, mean was used') % confirm
    
    elseif mostnans < tinygap
        fl = fillmissing(fl,'linear'); 
%         disp('ALL: TINY gap, linear was used') 
    
    end
    
    
    % fix gap after summer
    PostSummerGap = num2str(max(sum(ismissing(fl((sumrow:end),:)))));
    if PostSummerGap > 0
        gapnan = max(sum(ismissing(fl))); % Count NaN values in each column
        fl = fillmissing(fl,"movmean",hours(gapnan+1));
        % check for gaps
        gapnan = max(sum(ismissing(fl))); % check again
        count = 0; % keep count just for fun of the visualization
    
        while gapnan > 0 % run this until gapnan becomes zero
            fl = fillmissing(fl,"movmean",hours(gapnan+1)); % patch again
            gapnan = max(sum(ismissing(fl))); % check again
            count = count + 1; % count iterations
%             disp(['Iteration ' num2str(count) ', Current value of gapnan: ', num2str(gapnan)]);
        end
    
        if gapnan == 0 % final check if there are any NaNs left
%             disp('Post-summer patched; please check')
        end
    
    else
%         disp('No NaN gaps post-summer; please double-check')
    end

    nantest = num2str(max(sum(ismissing(fl(:,:)))));
    disp([nantest ' NaNs in ' num2str(yy)])

    assignin('base',char(nom),fl); % save separately
    
end

clear nanCount resultTable mostnans ans biggap count gapnan PostSummerGap...
    sumnans sumrow tinygap yy nom nantest fl maxhole avghole totalholes t1 t2...
    nan_lengths current_length porciento totalnans albb sumalbb i


%% SAVED AS:

load([mp 'drive_python/outputs/runoff_temp_feedback/' aws '/read_Promice_AWS/mets_2009_2022'])






%% 5) Plot check for albedos

% Double check where you save your figures
% FolderPath = [mp 'drive_python/figures/21. new albs oct11/' aws '/']; % location

% загрузи таблицу НаНов
load([mp 'drive_python/outputs/manuscript/albedo/met/' aws '/naninfo_table.mat'])

for yy = sy:2022
    % подготовь файл
    fl = ['met' num2str(yy)]; % get the met name
    fl = eval(fl); % get the actual timetable info  
    x = fl.Time; y1 = fl.albedo; 

    % Подготавливаем подписи на основе данных из таблицы
    annotations = cell(1,size(naninfo_table,1));
    for i = 1:size(naninfo_table, 1)
            row_label = naninfo_table.Properties.RowNames{i};
            row_value = naninfo_table{i, yy - (sy - 1)};
            annotations{i} = [row_label, ': ', num2str(row_value{1, 1})]; % Combine label and value
    end

    % рисуем
    figure(yy-(sy-1)); plot(x,y1); ylabel('albedo'); hold on;
    title([aws ' albedo (' num2str(yy) ')'])

    % добавь серый прямоугольник чтобы обозначить лето
    start_date = datetime(yy,6,1,0,0,0); end_date = datetime(yy,8,31,23,0,0);
    summer_indices = find(x >= start_date & x <= end_date);
    area([x(summer_indices);flip(x(summer_indices))],[ones(size(summer_indices));zeros(size(summer_indices))], ...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',1); % Set FaceAlpha to 1
    ylim([0.2 1]); hold off

    % Calculate the position for the annotations inside the graph area
    graphWidth = 0.8; % Adjust based on your graph dimensions
    graphHeight = 0.8; % Adjust based on your graph dimensions
    xPosition = 10 * graphWidth; % You can adjust the horizontal position
    yPosition = 0.61 * graphHeight; % You can adjust the vertical position

    % Add the title for summer stats
    title_text = 'Summer (shaded area) stats:';
    text(...
        xPosition, yPosition, title_text, ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'none', ...
        'EdgeColor', 'black', ...
        'FontSize', 10 ...
    );
    yPosition = yPosition - 0.05; % Adjust the vertical spacing

    % Add annotations as text boxes inside the graph area
    for i = 1:numel(annotations)
        text(...
            xPosition, yPosition, annotations{i}, ...
            'HorizontalAlignment', 'left', ...
            'BackgroundColor', 'none', ...
            'EdgeColor', 'none', ...
            'FontSize', 9 ...
        );
        yPosition = yPosition - 0.05; % Adjust the vertical spacing
    end



     % save figure
%     nom = [aws num2str(yy) '.png']; % save as this file name
%     FullPath = fullfile(FolderPath,nom);
%     print(FullPath, '-dpng', ['-r' num2str(300)]); % save with dpi 300
end

clear FolderPath FullPath x y1 yy nom fl i annotations row_value row_label...
    naninfo_table yPosition xPosition graphHeight graphWidth

% close all

% finished working on this section on oct 16, 2023, в 9:21




%% 6) Compare parameters of kanm to kanl ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% load and rename kanm met files
load([mp 'drive_python/outputs/runoff_temp_feedback/KANM/read_Promice_AWS/mets_2009_2022'])

for yy = 2009:2022
    met = ['met' num2str(yy)]; % get the met name
    met = eval(met); % get the actual timetable info
    nom = ['kanm' num2str(yy)]; % rename; have kanm instead of met
    assignin("base",char(nom),met) % save var into workspace
    varName = sprintf('met%d', yy); % id files that start with met
    clear(varName); % clear em
end

% kanl
load([mp 'drive_python/outputs/runoff_temp_feedback/KANL/read_Promice_AWS/mets_2009_2022'])
for yy = 2009:2022
    met = ['met' num2str(yy)]; % get the met name
    met = eval(met); % get the actual timetable info
    nom = ['kanl' num2str(yy)];
    assignin("base",char(nom),met)
    varName = sprintf('met%d', yy); % id files that start with met
    clear(varName); % clear em
end
clear yy nom varName met

% kanl
load([mp 'drive_python/outputs/runoff_temp_feedback/KANU/read_Promice_AWS/mets_2009_2022'])
for yy = 2010:2022
    met = ['met' num2str(yy)]; % get the met name
    met = eval(met); % get the actual timetable info
    nom = ['kanu' num2str(yy)];
    assignin("base",char(nom),met)
    varName = sprintf('met%d', yy); % id files that start with met
    clear(varName); % clear em
end
clear yy nom varName met

%% check for 2021

fm = kanm2021;
fl = kanl2021;
fu = kanu2021;
% Get the variable names of the timetables
kanlVars = fl.Properties.VariableNames;
kanmVars = fm.Properties.VariableNames;
kanuVars = fu.Properties.VariableNames;

x = fm.Time; % for plots in the loop
FolderPath = [mp 'drive_python/figures/15. kan_var_comparison_may29/2021_kan_all']; % saving location

% Loop over the variable names
for i = 1:numel(kanmVars)
    varName = kanmVars{i};
    
    % i) Compare the corresponding variable in kanl2021 and kanm2021
%     kanlData = fl.(varName);
    kanmData = fm.(varName);
    kanuData = fu.(varName);
    kanlData = fl.(varName);
    
    % ii) make plots to compare
    figure(i); % create new figure
    plot(x,kanlData,'DisplayName','kanl')
    hold on
    plot(x,kanmData,'DisplayName','kanm')
    plot(x,kanuData,'DisplayName','kanu','Color',"#7E2F8E")
    xlabel('Time'); ylabel(varName); legend
    hold off

    % iii) save figures
    nom = ['kan_um_' varName '.png']; % save as this file name
    FullPath = fullfile(FolderPath,nom);
    print(FullPath, '-dpng', ['-r' num2str(300)]); % save with dpi 300

end

close all
clear fm fl kanlVars kanmVars FolderPath x kanlData kanmData i varName nom...
    FullPath kanuData kanuVars






