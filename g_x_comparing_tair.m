% Comparing air temperatures (2m) for AWS and models
% Started on Oct 27, 2023, Maxim A Shapovalov





%% 0) Преамбула
clear; clc; % clear workspace and command window for neatness
% make future paths shorter with a shortcut (mp = max's path)
mp = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/'; 
% specify models and AWSs
aws = {'KANL' 'KANM' 'KANU'};
data = {'met' 'merra' 'mar' 'racmo'};
% data = {'met' 'mar'}; % for now (overwrite)

% % establish colors for the plotting
% colors = ({[0 0 0]; % met (black)
%     [0.9 0.5 0.1]; % merra (red-orange)
%     [0 0.4470 0.7410]; % mar (dark blue)
%     [0.4010 0.8450 0.9330]})'; % racmo (light blue)

% establish colors for the plotting
colors = ({[0 0 0]; % modis (black)
    [0.9 0.5 0.1]; % merra (red-orange)
    [0 0.4470 0.7410]; % mar (dark blue)
%     [0.4010 0.8450 0.9330]})'; % og racmo (light blue)
    [0.6350 0.0780 0.1840]})'; % new racmo (maroon)

% establish start and end years
sy = 2010; ey = 2021;


%% 1) Загрузка 

for a = 1:length(aws) % first, loop through weather stations
    met = aws{a}; % choose current MET station

    for b = 1:length(data) % second, loop through different models
        % 1) load data
        if b == 1 % if met
            load([mp 'drive_python/outputs/runoff_temp_feedback/' met '/read_Promice_AWS/mets_2009_2022'])
            % extract just tair
            for yr = sy:ey
                fl = eval(['met' num2str(yr)]); % get current file
                z = fl(:,1); % save just tair;
                nm = ['met' num2str(yr)]; % name of file
                assignin('base',nm,z); % create variable
                clear met2022 z nm fl
            end
            clear yr
        else % if merra, mar, or racmo
            load([mp 'drive_python/outputs/manuscript/tair/' data{b} '/' met '/' data{b} '_tair.mat'])
        end

    end
    % Now all files are loaded in. Begin Comparison!

    % 1) JJA scatter plots (met vs model)

    for model = 1:length(data)

        % prepare to extract JJA data
        temp = eval([data{model} num2str(sy)]);
        tiempo = temp.Time; % the time variable
        sstart  = datetime(year(tiempo(1)), 6, 1, 0, 0, 0); % june 1
        eend = datetime(year(tiempo(end)), 8, 31, 23, 0, 0); % aug 31
        rowsInPeriod = tiempo >= sstart & tiempo <= eend; % Find the rows that fall within the specified time period
        jja_rows = sum(rowsInPeriod == 1); % how many summer rows
        % Create arrays for storage
        dataCell = nan(jja_rows,length(sy:ey)); % all values will go here
        data_avg = nan(jja_rows,1); % the averages will go here

        for yr = sy:ey
            fl = eval([data{model} num2str(yr)]); % evaluate current file
            fljja = fl(rowsInPeriod,:); % extract just the summer data
            dataCell(:,(yr-(sy-1))) = fljja{:,1};
        end
        % get the averages
        for yr = 1:height(dataCell) 
            data_avg(yr,1) = mean(dataCell(yr,1));
        end

        % export to workspace the current data_avg file
        namee = [data{model} '_' met '_avg'];
        assignin('base',char(namee),data_avg);

    end

    % clear out the files
%     for yy = sy:ey
%         for varNames = data
%             fullname = [varNames{1},num2str(yy)];
%             clear(fullname)
%         end
%     end
    clear b data_avg dataCell eend fl fljja fullname jja_rows model...
        namee rowsInPeriod sstart tiempo varNames yr yy temp



   
    % 2) yearly line plots, individually vs met % % % % % % % % % % % % % % % 
%     for yr = sy:ey
%         metfl = eval(['met' num2str(yr)]); % get current met file
%         t = metfl.Time; % get variable for time for ease
%         for model = 2:4
%             otrofl = eval([data{model} num2str(yr)]); % get other file
%             figure; % create new figure
%             plot(t,metfl.tair,'DisplayName',met,'Color',colors{1}); hold on;
%             plot(t,otrofl.tair,'DisplayName',data{model},'Color',colors{model}); legend;
%             xlabel('Time'); ylabel('2-m Air Temperature [K]'); grid on;
%             title(['Comparing Air Temp of MET (' met ') and ' data{model} ' [' num2str(yr) ']'])
%             datetick('x','mmm','keepticks'); % now the year is not displayed (хорошо)!
% 
%             % save figure
% %             figFolderPath = [mp 'drive_python/figures/24. tair/1_YearlyVsMet/' met '/']; % location
% %             nom = [data{model} num2str(yr)]; % save as this file name
% %             FullPath = fullfile(figFolderPath,nom);
% %             print(FullPath, '-dpng', ['-r' num2str(300)]); % save with dpi 300
% 
%         end
%         clear model metfl t otrofl
%     end

    % 3) average plots (2010-2021)
    % i) Extract tair for all years and put together for each dataset
    for model = 1:length(data)
        dataCell = nan(height(mar2010),length(sy:ey)); % prep for storage
        data_avg = nan(height(mar2010),1); % the averages will go here
        % put together
        for yr = sy:ey
            fl = eval([data{model} num2str(yr)]);
            dataCell(:,(yr-(sy-1))) = fl{:,1};
        end
        % get the averages
        for yr = 1:height(dataCell) 
            data_avg(yr,1) = mean(dataCell(yr,1));
        end
        % plot
        figure(a); % fig for each aws
        plot(mar2010.Time,data_avg,'DisplayName',data{model},'Color',colors{model})
%         title(['Comparing Average Air Temperatures for ' met ' (2010-2021)'])
        xlabel('Time'); ylabel('2-m Air Temperature [K]')
        datetick('x','mmm','keepticks'); % now the year is not displayed (хорошо)!
        hold on; grid on; ylim([230 290])
        % specific aws legend
        if a==1 % if kanl
            legend('KAN-L','MERRA-2','MAR','RACMO');
        elseif a==2 % if kanm
            legend('KAN-M','MERRA-2','MAR','RACMO');
        else % if kanu
            legend('KAN-U','MERRA-2','MAR','RACMO');
        end
        % assign as distinct name in workspace
        nom = [data{model} '_avg']; assignin('base',nom,data_avg);
        % save plot
%         figFolderPath = [mp 'drive_python/figures/24. tair/2_AvgYears/']; % location
%         nom = ['avg_Air_Temps_' met]; % save as this file name
%         FullPath = fullfile(figFolderPath,nom);
%         print(FullPath, '-dpng', ['-r' num2str(300)]); % save with dpi 300
    end

    % Calculate RMSEs
    mar1 = (met_avg - mar_avg).^2; rmse_mar = sqrt(mean(mar1));
    disp([data{3} 'rmse for ' met ' is ' num2str(rmse_mar)])
    merra1 = (met_avg - merra_avg).^2; rmse_merra = sqrt(mean(merra1)); 
    disp([data{2} 'rmse for ' met ' is ' num2str(rmse_merra)])
    racmo1 = (met_avg - racmo_avg).^2; rmse_racmo = sqrt(mean(racmo1)); 
    disp([data{4} 'rmse for ' met ' is ' num2str(rmse_racmo)])


end

%% Opt 1 continued: scatter plots

% Create a tiled layout
tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'normal');
% title for all 9 panels
% sgtitle('Correlating average 2-m air temperatures [K] of METs to other models (2010-2021)','FontWeight','Bold','FontSize',14);
set(gcf,'Position',[50, 75, 600, 500]) % figure size
ScatterData = {'merra' 'mar' 'racmo'};


for b = 1:length(ScatterData)

    for a=1:length(aws)

        nexttile(a + (b-1)*3); % start plotting things on the tiles(subplots)

        met = aws{a};

        xnom = [data{1} '_' aws{a} '_avg']; % select the x
        ynom = [ScatterData{b} '_' aws{a} '_avg']; % and the y (just names)
        x = eval(xnom); y = eval(ynom); % get the specific x and y

        scatter(x,y,7,colors{b+1},'filled'); hold on; % make scatter plot
        % Create and add linear fit
        p = polyfit(x, y, 1); y_fit = polyval(p, x);
        plot(x, y_fit, 'r-', 'DisplayName', 'Linear Fit');
        % Caclulate and add RMSE
        rmse = sqrt(mean((y_fit - y).^2)); 
        formatted_rmse = sprintf('%.2f',rmse); % control how many decimal points
        text(0.6, 0.05, ['RMSE: ' num2str(formatted_rmse)], 'Units', 'normalized');
        xlim([260 290]); ylim([260 290]); grid on;



        if b==1 && a==1
            ylabel('MERRA-2'); title('KAN-L','FontWeight','normal')
            xticklabels([]);
        elseif b==1 && a==2
            title('KAN-M','FontWeight','normal'); xticklabels([]); yticklabels([]);
        elseif b==1 && a==3
            title('KAN-U','FontWeight','normal'); xticklabels([]); yticklabels([]);
        elseif b==2 && a==1
            ylabel('MAR (Temperature [K])'); xticklabels([]);
        elseif b==2 && a==2
            xticklabels([]); yticklabels([]);
        elseif b==2 && a==3
            xticklabels([]); yticklabels([]);
        elseif b==3 && a==1
            ylabel('RACMO')
        elseif b==3 && a==2
            xlabel('AWS (Temperature [K])'); yticklabels([]);
        else
            yticklabels([]);
        end


        hold off;

    end


end








