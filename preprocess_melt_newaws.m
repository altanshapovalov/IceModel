% Pre-process data for scripts that create manuscript figures (newer aws)
% Purpose: put data after icemodel in here, save the output, and put that output
% into scripts fig2_ fig3_ etc that create figs for manuscript
% Created march 19, 2024; Maxim Altan-Lu Shapovalov

%% Preamble

clear; clc; % clear workspace and command window for neatness
% make future paths shorter with a shortcut (mp = max's path)
mp = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/'; 
mpn = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel_new/';
data = {'modis' 'merra' 'mar' 'racmo'}; % albedo sources
% data = {'modis'}; % just for getting mcd43a3 

% % specify models and AWSs
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
awsz = {'kanl','kanm','kanu','kpcl','kpcu','nukl','nuku','qasl','qasu','scol',...
    'scou', 'thul', 'thuu', 'upel', 'upeu'};

% save figs here
% FolderPath = [mp 'drive_python/figures/29. new aws/15. avg melts/'];
% FolderPath = [mp 'drive_python/figures/29. new aws/34. mcd43a3 single pixel/4. avg melts/'];
% save avg albedo data here
% outloc = [mp 'drive_python/outputs/manuscript/melt/new_aws/mcd43a3/single pixel/']; 
outloc = [mp 'drive_python/outputs/manuscript/melt/newaws_w_iqr/']; 




%% 1) Organize all melts for given AWS in a table

for xx = 1:length(awsz)
    aws = awsz{xx};
    yrz = specifyYears(aws); % use function to select specific aws years

    for b = 1:length(data) % loop through different models
        % i) prepare empty array to store albedo of all years - - - - - - - - - - - 
        melts(8760,length(yrz)) = NaN; % 8760 bc non leap year # of hours
        melts(:,:) = NaN; % ensures that full of NaNs and not zeros
        % ii) load data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        for y1 = 1:length(yrz) % loop through each available year
            % a) load melt (specific year, aws, and data source)
            if b==1 % if modis (mcd43a3)
                load([mpn 'output/mcd43a3/single pixel/ice1/' aws '_' data{b} '_' num2str(yrz{y1}) '.mat'])
            else
                load([mpn 'output/ice1/' aws '_' data{b} '_' num2str(yrz{y1}) '.mat'])
            end
            % b) save into albss
            melts(:,y1) = ice1.melt;
        end
        % iii) get average (median) for each row/hour across all years - - - - - - -
%         avgg(height(melts),3) = NaN; avgg(:,:) = NaN; % prep array for storage
        avgg(height(melts),1) = NaN; avgg(:,:) = NaN;
        iqr75 = avgg; iqr25 = avgg; % prep the iqr vars too
        iqrr = avgg; % and prep this
        for c = 1:height(avgg)
            avgg(c) = median(melts(c,:));
            iqr75(c) = prctile(melts(c,:),75);
            iqr25(c) = prctile(melts(c,:),25);
%             avgg(c,:) = quantile(melts(c,:),[.25 .50 .75]);
            tv1 = iqr75(c) - avgg(c); % see what the diff is btwn 75 and 50%
            tv2 = avgg(c) - iqr25(c);
            tv3 = (tv1 + tv2) / 2; % get avg of the two for my single iqr+-
            iqrr(c) = tv3;
        end
        avg1.Time = ice1.Time; % put in a structure
%         avg1.melt25 = avgg(:,1); 
%         avg1.melt50 = avgg(:,2); 
%         avg1.melt75 = avgg(:,3); 
        avg1.melt = avgg;
        avg1.iqr = iqrr;
        avg = table2timetable(struct2table(avg1)); % convert to timetable
    %     avg = retime(avg,'daily',@median); % convert to daily
        % iv) save file - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        newname = [aws '_' data{b}]; % prep unique name
        nnal = [outloc newname]; % new name and location
        save(nnal,'avg'); % save data in new location
        % v) create and save plot(s) - - - - - - - - - - - - - - - - - - - - - - - -
%         figure('visible', 'off'); % create figure without displaying it
%         for t = 1:length(yrz) % plot each year
%             plot(ice1.Time,melts(:,t),'Color','k'); hold on;
%         end
%         plot(ice1.Time,avgg(:,2),'Color','r','LineWidth',2.5) % plot the avg over it
%         title(['Average ' data{b} ' melt at ' aws]); % unique title
%         xlabel('Time'); ylabel('Melt [m]'); 
%         datetick('x','mm','keepticks'); % don't print the year (since it's average)
%         nom = [aws '_' data{b} '.png']; % specify name
%         FullPath = fullfile(FolderPath,nom); % establish full path
%         print(FullPath,'-dpng',['-r' num2str(300)]); % save with 300 dpi

        disp([aws ' finished'])
        disp(' ')
    end
    clear avg1 melts y1 c b newname nnal
end
disp('Финита ля комедия!')




