% compare single pixel mcd43a3 albedo to 5 and 10 km resolutions
% Created by Maxim Altan-Lu Shapovalov on May 21, 2024

%% Преамбула
clear; clc
mp = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/'; % mp = max's path

FolderPath = [mp 'drive_python/figures/29. new aws/34. mcd43a3 single pixel/6. compare resolutions/']; % save figs here

opt = 'modis'; % comparing modis albedos

% CHANGE THIS FOR DIFFERENT WEATHER STATION
% aws = 'kanl';
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
% aws = {'kanl','kanm','kanu','kpcl','kpcu','nukl','nuku','qasl','qasu','scol',...
%     'scou', 'thul', 'thuu', 'upel', 'upeu'};

% select aws (ones that are far enough from the GrIS margin)
aws = {'kanm','kanu','kpcu','qasu','upeu'}; % maybe also nuku

% Что мы сравниваем с единичным пикселем?
% sravneniye = '5km'; % mcd43a3 5km
sravneniye = '10km'; % mcd43a3 10km

%% 1) Load albs

for a = 1:length(aws) % loop thru weather stations
    % i) specify station and years
    met = aws{a}; % choose current aws
    [yrz,elev,latt,lonn] = specifyYears(met); % get years

    for b = 1:length(yrz) % loop thru the years

        % ii) load albs
        load([mp 'input/userdata/mcd43a3/single pixel/' opt '_' met '_' num2str(yrz{b}) '.mat'])
        sp = Data.albedo; % save the single pixel albedo
        clear Data % на всякий пожарный
        % choose options:
        if strcmp(sravneniye,'5km')
            load([mp 'input/userdata/mcd43a3/' opt '_5km_' met '_' num2str(yrz{b}) '.mat'])
            srav_alb = Data.albedo;
        elseif strcmp(sravneniye,'10km')
            load([mp 'input/userdata/mcd43a3/' opt '_10km_' met '_' num2str(yrz{b}) '.mat'])
            srav_alb = Data.albedo;
        end
        T = Data.Time; % save the time
        clear Data % just to be sure
        % Get just summer info
        row1 = 3625; % june 1, 00:00
        row2 = 5832; % aug 31, 23:00
        spsum = sp(row1:row2,1); % just summer albedo (mod)
        sravsum = srav_alb(row1:row2,1); % mcd
        % Apr1-oct1
        ap1 = 2161; % apr1 00:00
        oc1 = 6553; % oct1 00:00
        sp_aproct = sp(ap1:oc1,1);
        srav_aproct = srav_alb(ap1:oc1,1);

        % iii) compare
        % Create a scatter plot
        figure('Visible','off') % don't print every figure
        plot(T,sp,'DisplayName','Single Pixel'); hold on; 
        if strcmp(sravneniye,'5km')
            plot(T,srav_alb,'DisplayName','5 km');
        elseif strcmp(sravneniye,'10km')
            plot(T,srav_alb,'DisplayName','10 km');
        end
        xlabel('Time');
        ylabel('Albedo [-]');
%         title(['5 vs 10 km resolutoions of MCD43A3 at ' met ' (' num2str(yrz{b}) ')']);
        ylim([0.3 1.0])
        legend
        
        % Perform linear regression + + +
        p = polyfit(sp, srav_alb, 1);
        yfit = polyval(p, sp);
        % Calculate the R-squared value
        yresid = srav_alb - yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(srav_alb)-1) * var(srav_alb);
        rsqf = 1 - SSresid/SStotal;
        % same, but for summer + + + 
        p = polyfit(spsum, sravsum, 1); yfit = polyval(p, spsum);
        % Calculate the R-squared value
        yresid = sravsum - yfit; SSresid = sum(yresid.^2);
        SStotal = (length(sravsum)-1) * var(sravsum); rsqs = 1 - SSresid/SStotal;
        % same, but for apr-oct + + +
        p = polyfit(sp_aproct, srav_aproct, 1); yfit = polyval(p, sp_aproct);
        % Calculate the R-squared value
        yresid = srav_aproct - yfit; SSresid = sum(yresid.^2);
        SStotal = (length(srav_aproct)-1) * var(srav_aproct); rsqao = 1 - SSresid/SStotal;

        % Display R-squared value on the plot
%         text(5, 0.97, ['Full Year R^2 = ', num2str(rsqf)], 'FontSize', 12, 'Color', 'k');
        text(5, 0.97, ['Apr-Oct R^2 = ', num2str(rsqao)], 'FontSize', 12, 'Color', 'k');
        text(5, 0.93, ['Summer R^2 = ', num2str(rsqs)], 'FontSize', 12, 'Color', 'k');
        hold off;

        % save fig
        nom = [met '_' sravneniye '_' num2str(yrz{b}) '.png']; % specify name
        FullPath = fullfile(FolderPath,nom); % establish full path
        print(FullPath,'-dpng',['-r' num2str(300)]); % save with 300 dpi
        
    end
    disp(['Finished saving ' met])
    disp(' ')

end
disp('Touts sont finis!')