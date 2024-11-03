% Preparing melt data for figures
% Created on Jan 31, 2024; Maxim Altan-Lu Shapovalov

% Purpose: icemodel 1.0 provided ice1 timetable (that has melt) for each year
% (2010-2021) and dataset (mar, merra, modis, racmo). Must extract melt,
% organize it in a table, confirm correct numbers were saved, obtain averages
% over study period (avg of Jan 1, avg of Jan 2, etc). Must do this for each
% dataset at each AWS site (KANL,M,U). Files will then be read by scripts titled
% "fig_" in scripts>manuscript.


% ********
% IMPORTANT: code was run on jan31 & feb1, 2024 (AND on feb 8, 2024 for
% quantiles). If run again, files can be overwritten, so be aware of what you 
% do with this script.
% ********

clear; clc % clear workspace and command window
% parameters
sitename = 'kanu';        % options: 'kanm' 'kanl' 'kanu'
userdata = 'mar';       % options: 'modis','racmo','merra','mar'
simyears = 2010:2021;     % study period (simulation years)
forcings = sitename;      % the forcing will always match the sitename
mp = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/'; % mp = max's path
mpn = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel_new/'; % mpn = max's path (new icemodel)


%% 0) Extract melt and store

% a) clarifications and loading of file
outp = [mpn 'output/ice1/']; % specify path to output files (with ice1)
i = simyears(1); % extract a timetable just to use its time
filename = [sitename '_' userdata '_' num2str(i)]; % specify which file
load([outp filename '.mat']) % load the specified file
tiempo = ice1.Time; % will be used in the future

% b) create table in which values will be stored
numRows = height(ice1); numCols = numel(simyears); % dimensions of storage table
dataArray = zeros(numRows,numCols); % create table for storage
disp(['Performance for ' userdata ' at ' sitename]); % be aware

% c) place all data into a single array for organization
for k = simyears
    filename = [sitename '_' userdata '_' num2str(k)]; % specify which file
    load([outp filename '.mat']) % load the specified file
    dataArray(:,k-(simyears(1)-1)) = ice1.melt;
    disp([num2str(k) ' saved into column ' num2str(k-(simyears(1)-1))])
    % also, create individual timetables for fig 5
    uniquename = [userdata num2str(k)]; % prepare name that fig5 code can read
    themelt = ice1.melt; % save melt into a variable
    thedata.Time = ice1.Time; thedata.melt = themelt; % put into structure
    ohyeah = table2timetable(struct2table(thedata)); % convert to timetable
    assignin('base',uniquename,ohyeah); % create individual array
    % save all 12 timetables manually; saved in output>all_melts
end
clear ohyeah thedata themelt
% quick test to see if the data matches
% corr(dataArray(:,12),mar2021.melt) % change the second file name accordingly
% % export as table for fig 5 (obsolete); went into allmelts>in_single_array
% dataTable = array2table(dataArray);
% dataTable.Properties.VariableNames = cellstr(string(simyears));
% nom = [sitename '_' userdata '_allmelts'];
% path2 = [mpn 'output/allmelts/'];
% exppath = [path2 nom '.mat'];
% save(exppath,'dataTable')


%% d) test if values match (i.e., if correctly transferred data)
% load select year
yr = 2017; nm = [sitename '_' userdata '_' num2str(yr)]; load([outp nm '.mat'])
wellfit = corr(ice1.melt,dataArray(:,8)); disp(['Correlation: ' num2str(wellfit)])

% % e) get average by day across study years
avgg = zeros(numRows,1); % prepare empty array
for c = 1:numRows
    avgg(c) = mean(dataArray(c,:));
end

% e2) get 25th and 75th quantiles
% e) get average by day across study years
qtls = zeros(numRows,3); % prepare empty array (melt 25th quantile)
for c = 1:numRows
    qtls(c,:) = quantile(dataArray(c,:),[.25 .50 .75]);
end

% PLOT TO COMPARE MEAN TO MEDIAN AND .25 + .75 QUANTILES
% lww = 2; % line width
% xx = tiempo; % time from timetable
% % plot(xx,avgg(:,:),'Color',"#77AC30",'LineWidth',lww,'DisplayName','Mean'); hold on;
% % plot(xx,qtls(:,2),'Color',"#7E2F8E",'LineWidth',lww,'DisplayName','Median'); legend; 
% % xlabel('Time'); ylabel('Melt [m]'); 
% % title(['Mean vs Median for ' userdata ' at ' sitename ' (2010-2021)']); hold off
% % Plot 2 – plot all quantiles and mean together
% plot(xx,qtls(:,1),'Color','b','LineWidth',lww,'DisplayName','25th'); hold on
% plot(xx,qtls(:,2),'Color',"#7E2F8E",'LineWidth',lww,'DisplayName','Median')
% plot(xx,qtls(:,3),'Color','r','LineWidth',lww,'DisplayName','75th'); 
% plot(xx,avgg(:,:),'Color',"#77AC30",'LineWidth',lww,'DisplayName','Mean');
% legend('Location','west'); ylabel('Melt [m]'); xlabel('Time');
% title(['Quartiles for ' userdata ' at ' sitename ' (2010-2021)']); hold off
% % set(gcf,'Position',[600, 500, 1000, 350]) % figure size
% datetick('x','mmm','keepticks'); % now the year is not displayed (хорошо)!
% 
% % save figure
% savepath1 = [mp 'drive_python/figures/28. paper1_v3/melt quantiles feb8/'];
% figname = [sitename '_' userdata '_melt_quantiles'];
% savepath2 = [savepath1 figname '.png'];
% print(savepath2, '-dpng',['-r' num2str(300)]); % save with dpi 300


%%
% f) assign into a timetable with proper name
% finalfile.Time = tiempo; finalfile.melt = avgg; % assign into a structure
finalfile.Time = tiempo; 
finalfile.twentyfifth = qtls(:,1); % assign into a structure
finalfile.median = qtls(:,2); % assign into a structure
finalfile.seventyfifth = qtls(:,3); % assign into a structure

final = table2timetable(struct2table(finalfile)); % convert to timetable
% newname = [sitename '_' userdata '_avgmelt']; % for mean
newname = [sitename '_' userdata '_melt_quantiles']; % for quantiles
assignin('base',newname,final);

% g) export
path1 = [mpn 'output/avgmelts/'];
exportpath = [path1 newname '.mat'];
save(exportpath,'final')


% % in case of leap year:
% if height(ice1) == 8784 % if leap year
%     feb29 = month(ice1.Time) == 2 & day(ice1.Time) == 29; %id feb29
%     ice1 = ice1(~feb29,:); % update ice1, now w/ no leap day data
%     disp([num2strz(i) ' was corrected for leap data'])
% elseif height(ice1) ~= 8784 && height(ice1) ~= 8760
%     error('ice1.Time should be either 8760 or 8784; double check')
% end

