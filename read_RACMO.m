% RACMO data analysis (albedo & swsd)
% started on March 20, 2023

% RACMO2.3p2 at 5.5 km:
% Daily surface albedo (alb, -) and shortwave downward radiation (swsd, W/m2) for 1 January 2000 - 31 December 2022.
% Mask file including: lon/lat, PROMICE mask (4 = GrIS, 2 and 3 = GICs), grid-cell area and topography from the GIMP DEM upscaled to 5.5 km.
% Forcing: The data set is based on ERA5 (2000-2022) reanalyses.
% Reference: The data set is an extension of Noël et al. (2019).

% organizing and modifying on June 23, 2023

% Contents:
% ALBEDO ============
    % 1) read nc file
    % 2) extract data



%%
clear
mp = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/'; % mp = max's path

% CHANGE THIS FOR DIFFERENT WEATHER STATION"
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

% specify start and end years
sy = 2009;
ey = 2022;

% aws location
if strcmp(aws,'kanm') % if KANM
    coord_target = [-48.8329, 67.0668]; % get the KANM lon and lat
elseif strcmp(aws,'kanl') % if KANL
    coord_target = [-49.9482, 67.0957]; % kanl location
elseif strcmp(aws,'KANU') % if KANU
    coord_target = [-47.0255, 67.0007];
elseif strcmp(aws,'kpcl') % if KANU
    coord_target = [-24.0829, 79.9109];
elseif strcmp(aws,'kpcu') % if KPC_U
    coord_target = [-25.1666, 79.8345];
elseif strcmp(aws,'nukl') % if NUK_L
    coord_target = [-49.5337, 64.4829];
elseif strcmp(aws,'nuku') % if NUK_U
    coord_target = [-49.2714, 64.5063];
elseif strcmp(aws,'qasl') % if QAS_L
    coord_target = [-46.852, 61.03];
elseif strcmp(aws,'qasu') % if QAS_U
    coord_target = [-46.8195, 61.1752];
elseif strcmp(aws,'scol') % if SCO_L
    coord_target = [-26.8176, 72.2258];
elseif strcmp(aws,'scou') % if SCO_U
    coord_target = [-27.2578, 72.3942];
elseif strcmp(aws,'thul') % if THU_L
    coord_target = [-68.2667, 76.3998];
elseif strcmp(aws,'thuu') % if THU_U
    coord_target = [-68.1483, 76.4193];
elseif strcmp(aws,'upel') % if UPE_L
    coord_target = [-54.2955, 72.8932];
elseif strcmp(aws,'upeu') % if UPE_U
    coord_target = [-53.5748, 72.8882];
else
    error([aws ' is NOT a valid weather station on the K-transect' ...
        '(options: KANL, KANM, KANU, KPC_L, KPC_U, NUK_L, NUK_U, QAS_L, ' ...
        'QAS_M, QAS_U, SCO_L, SCO_U, THU_L, THU_U, UPE_L, UPE_U)']);
end




%% =================== ALBEDO ==================================================
for yr = sy:ey

    % 1) read nc file ----------------------------------------------------------
    % get path in which file is located
    pth = [mp 'drive_python/input/RACMO/Daily_5.5km/alb/alb.' num2str(yr) '.FGRN055_BN_RACMO2.3p2_ERA5_3h_FGRN055.DD.nc'];
    
    ncfile = dir(pth); % get basic info
%     ncdisp(ncfile.name); % to see ALL info of .nc file
    % INFO: lon, lat (both 438x566),alb(438x566x1x366)
    fnom = ncfile.name; % get file name
    disp(['Analyzing ' num2str(yr) ' for ' aws])


    % 2) extract data ----------------------------------------------------------
    % some variables need only be extracted once, do them first
    lon    = ncread(fnom,'lon');
    lat    = ncread(fnom,'lat');
    swsd   = ncread(fnom,'alb');
    sizee  = size(swsd); % size of the nc file
    
    % ------------Insung coming to the rescue: -----------------
    rwz = sizee(1); % # of rows
    clmz = sizee(2); % # of columns. Shape is rwz x clmz
    lon_flat     = reshape(lon, [rwz*clmz,1]); % flatten lon array
    lat_flat     = reshape(lat, [rwz*clmz,1]); % flatten lat array
    lonlat_flat  = [lon_flat, lat_flat];      % bring them together

    
    % see which lon/lat pair is closest to my given coordinates:
    euclidean_from_target = sqrt((lonlat_flat(:,1) - coord_target(1)).^2 +...
        (lonlat_flat(:,2) - coord_target(2)).^2);
    
    % id the row in which we have those closest coordinates:
    indx_min_euc = find(euclidean_from_target == min(euclidean_from_target));
    
    % see what are the coord and make sure that they're close enough:
    targ_coord = lonlat_flat(indx_min_euc,:); 
    
    indx_col_target = ceil(indx_min_euc / rwz); % find col that has these coords
    indx_row_target = indx_min_euc - (rwz * (indx_col_target - 1)); % find row
    addr_in_template = [indx_row_target, indx_col_target]; % row col in LAT LON
    racmo_alb = swsd(indx_row_target, indx_col_target, :, :); % get alb from our coords!
    rracmo(sizee(4),1) = NaN; % create empty array
    rracmo(:,1)   = racmo_alb(:,:,:,:); % fill it up!
    
    racmo_alb_daily = rracmo; % rename for better comprehension
    
    clear rwz clmz lonlat_flat lon_flat lat_flat euclidean_from_target...
        indx_min_euc indx_col_target indx_row_target addr_in_template racmo_alb...
        rracmo fnom lat lon swsd pth ncfile sizee % clear vars that won't be used later
    


   
    % 3) expansion -------------------------------------------------------------
    % a) turn it from daily to hourly 
    Nfiles = length(racmo_alb_daily); % how many days (leap or not)
    hrs    = 24; % hours in a day
    rac_alb((Nfiles*hrs),1) = NaN; % prep hourly array
    rac_alb(1:(length(rac_alb)),1) = NaN; % turn zeros to NaN
    
    % b) fill every day's 10:00 (am) row 
    for i = 1:Nfiles % for each day
    % NOTE: time will start w/ Jan1 00:00, so to put for 10:00, use row 11
        mar_row = 11; % data collected at 10:30 so will put as 10am
        indx      = mar_row + (hrs*(i-1)); % which row
        rac_alb(indx,1) = racmo_alb_daily(i,1); % insert value
    end
    clear Nfiles hrs i mar_row indx racmo_alb_daily
    
    % c) add time
    T = (datetime(yr,1,1):hours(1):datetime(yr,12,31)+hours(23))'; % hourly
    R.Time = T; R.albedo = rac_alb; % create structure
    racmo_1h = table2timetable(struct2table(R));
    clear rac_alb T R
    
    % d) if feb 29, delete it
    feb29 = month(racmo_1h.Time) == 2 & day(racmo_1h.Time) == 29;
    racmo_1h = racmo_1h(~feb29,:); % update mod, now w/ no leap days
    clear feb29
    
    % e) fill gaps
    racmo_1h.albedo = fillmissing(racmo_1h.albedo,'linear');
    
    % f) save albedo for year
    nom = (['racmo' num2str(yr)]); % specify file name
    assignin('base',char(nom),racmo_1h); % save separately

    % g) keep track of progress
%     disp([num2str(yr) ' done'])

end

clear coord_target racmo_1h

% realized that 2000-2008 won't be used. oh well
% SAVED AS:
% load([mp 'drive_python/outputs/albedo/RACMO/' aws '/racmo_albs_2000_2022.mat'])


%% Save the albedo so that IceModel v.1 can read it (added March 7th, 2024)
% Thus, rename each file to "Data", keep as timetable, save with naming
% convention of model_aws_year.mat (e.g., racmo_kpcl_2013.mat).
% Model options: modis, merra, mar, racmo.
% Save to icemodel>input>userdata

% Before that, save figures of albedos just to visually inspect them.
% In figures>29>6
for yyyy=sy:ey
    disp(['Saving ' num2str(yyyy) ' for ' aws])
    fl = eval(['racmo' num2str(yyyy)]); % get current file
    metVars = fl.Properties.VariableNames; % get var names
    tiempo = fl.Time; % the time
    FolderPath = [mp 'drive_python/figures/29. new aws/8. racmo albs/']; % save figs here
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
    nom = ['racmo' num2str(i)]; % get name of file
    Data = eval(nom); % activate select file
    outloc = [mp 'input/userdata/']; % output location (save to this new place)
    newname = ['racmo_' aws '_' num2str(i) '.mat'];
    nnal = [outloc newname]; % new name and location (for met file)
    save(nnal,'Data'); % save in the new location
end
clear i nom outloc newname nnal Data







% %% Plot Comparisons of MET, MODIS, and RACMO
% 
% % load racmo and mets for all years
% load([mp 'drive_python/outputs/albedo/RACMO/' aws '/racmo_albs_2000_2022.mat'])
% load([mp 'drive_python/outputs/albedo/MODIS/' aws '/modis_albs_2009_2022.mat'])
% load([mp 'drive_python/outputs/runoff_temp_feedback/' aws '/read_Promice_AWS/mets_2009_2022'])
% 
% 
% 
% % what would you like to compare, sir? [PLS CHOOSE ONE] * * * * * * * * * * * *
% opt1 = "met and racmo";
% % opt2 = "modis and racmo";
% % opt3 = "met, modis, and racmo";
% 
% 
% 
% % saving location
% if exist('opt1','var') % if option 1 var is present in workspace
%     FolderPath = [mp 'drive_python/figures/13. aws_albs_09_22/' aws '/3_MetRacmo']; % location
% elseif exist('opt2','var')
%     FolderPath = [mp 'drive_python/figures/13. aws_albs_09_22/' aws '/4_ModisRacmo']; % location
% elseif exist('opt3','var')
%     FolderPath = [mp 'drive_python/figures/13. aws_albs_09_22/' aws '/5_MetModisRacmo']; % location
% end
% 
% % establish the start year (sy)
% if strcmp(aws,'KANM') || strcmp(aws,'KANL') % KANL and KANM start with 2009
%     sy = 2009;
% elseif strcmp(aws,'KANU') % KANU starts with 2010
%     sy = 2010;
% end
% 
% for yy = sy:2022
% 
%     % a) get file
%     met = ['met' num2str(yy)]; % get the met name
%     met = eval(met); % get the actual timetable info
%     % same with modis and racmo now
%     modis = ['modis' num2str(yy)]; modis = eval(modis);
%     racmo = ['racmo' num2str(yy)]; racmo = eval(racmo);
% 
%     % b) set parameters
%     lw = 3; % linewidth
%     metcvet = 'k'; % noir pour met
%     modiscvet = "#EDB120"; % jaune pour modis
%     racmocvet = "#77AC30"; % vert pour racmo
%     ftsz = 17.5; % fontsize
%     ftnm = 'Damascus'; % fontname; alt Times, Arial
% 
%     % c) plot
%     if exist('opt1','var') % if option 1 var is present in workspace
%         plot(met.Time,met.albedo,'LineWidth',lw,'Color',metcvet,'DisplayName',['MET (' aws ')']); hold on;
%         plot(met.Time,racmo.albedo,'LineWidth',lw,'Color',racmocvet,'DisplayName','RACMO');
%     elseif exist('opt2','var')
%         plot(met.Time,modis.albedo,'LineWidth',lw,'Color',modiscvet,'DisplayName','MODIS'); hold on;
%         plot(met.Time,racmo.albedo,'LineWidth',lw,'Color',racmocvet,'DisplayName','RACMO');
%     elseif exist('opt3','var')
%         plot(met.Time,met.albedo,'LineWidth',lw,'Color',metcvet,'DisplayName',['MET (' aws ')']); hold on;
%         plot(met.Time,modis.albedo,'LineWidth',lw,'Color',modiscvet,'DisplayName','MODIS');
%         plot(met.Time,racmo.albedo,'LineWidth',lw,'Color',racmocvet,'DisplayName','RACMO');
%     end
% 
% 
%     % d) post-plot neatness // detail refinement
%     set(gcf,'Position',[50, 75, 1340, 725]) % figure size
%     % set title
%         if exist('opt1','var') % if option 1 var is present in workspace
%             title(['Comparing MET and RACMO albedo (' num2str(yy) ')'],'FontSize',ftsz+12,...
%             'fontname',ftnm); % title
%         elseif exist('opt2','var')
%             title(['Comparing MODIS and RACMO albedo at ' aws ' site (' num2str(yy) ')'],'FontSize',ftsz+12,...
%             'fontname',ftnm); % title
%         elseif exist('opt3','var')
%             title(['Comparing albedos of MET, MODIS, and RACMO (' num2str(yy) ')'],'FontSize',ftsz+12,...
%             'fontname',ftnm); % title
%         end
%     xlabel('Time','FontSize',ftsz,'fontname',ftnm); % x-axis name
%     ylabel('Albedo','FontSize',ftsz,'fontname',ftnm); % y-axis nom
%     lgd = legend; lgd.FontSize = 20; lgd.Location = 'east'; % legend
%     ax = gca; % control axes
%     ax.TitleFontSizeMultiplier = 1.5; % make the title larger
%     ax.LineWidth = 1.5;
%     set(gca,'fontsize',14); % enlargen size of axis numbers
%     ylim([0 1])
% 
%     % f) save figure
%     if exist('opt1','var') % if option 1 var is present in workspace
%         nom = ['MetRacmoAlb' num2str(yy) '.png']; % save as this file name
%     elseif exist('opt2','var')
%         nom = ['ModisRacmoAlb' num2str(yy) '.png']; % save as this file name
%     elseif exist('opt3','var')
%         nom = ['MetModisRacmoAlb' num2str(yy) '.png']; % save as this file name
%     end
%     FullPath = fullfile(FolderPath,nom);
%     print(FullPath, '-dpng', ['-r' num2str(300)]); % save with dpi 300
% 
%     hold off
% 
% end
% 
% close all
% 
% clear ans ax FolderPath ftnm ftsz FullPath lgd lower_limit lw metcvet minn...
%     modiscvet nom yy modis met racmocvet
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% 5) plot
% % y1 = racmo.albedo;
% % x = racmo.Time;
% % y = met.albedo; % import mar and met manually
% % y2 = mar.albedo;
% % plot(x,y,'LineWidth',2.5,'DisplayName','KAN-M')
% % hold on
% % plot(x,y1,'LineWidth',2.5,'DisplayName','RACMO')
% % plot(x,y2,'LineWidth',2.5,'DisplayName','MAR','Color','k')
% % 
% % 
% % set(gcf,'Position',[200, 200, 1500, 1000]) % controls figure size
% % fntsz = 15; % font size
% % ax = gca;  %get variable to control axes
% % fnt = 'Helvetica'; % font type; e.g., times, arial, helvetica, damascus
% % title('Comparison of KAN-M AWS and RACMO albedo for 2016')
% % ylabel('albedo','fontsize',fntsz,'fontname',fnt); 
% % xlabel('Time','FontSize',fntsz,'fontname',fnt);
% %     set(gca,'fontsize',fntsz); % enlargen size of axis numbers
% %     ax.TitleFontSizeMultiplier = 1.5; % make title larger
% %     ax.LineWidth = 1.5; grid on;
% % lgd = legend; lgd.FontSize = 16; lgd.Location = 'east';
% 
% 
% 
% 
% %%
% clear addr_in_template ans clmz coord_target euclidean_from_target indx_col_target...
%     indx_min_euc indx_row_target lat lat_flat lon lon_flat lonlat_flat pth rwz...
%     swsd racmo_alb_daily racmo_alb rracmo timez yr RA
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% =================== SWSD ===============================================
% 
% 
% %% 1) read nc file
% % pth = 'drive_python/input/RACMO/Daily_5.5km/alb/*.nc';
% pth = 'drive_python/input/RACMO/Daily_5.5km/swsd/swsd.2016.FGRN055_BN_RACMO2.3p2_ERA5_3h_FGRN055.DD.nc';
% 
% 
% ncfile = dir(pth); % get basic info
% ncdisp(ncfile.name); % to see ALL info of .nc file
% % INFO: lon, lat (both 438x566),alb(438x566x1x366)
% fnom = ncfile.name; % get file name
% 
% 
% 
% 
% %% 2) extract vars
% % some variables need only be extracted once, do them first
% lon    = ncread(fnom,'lon');
% lat    = ncread(fnom,'lat');
% swsd    = ncread(fnom,'swsd');
% 
% % ------------Insung coming to the rescue: -----------------
% rwz = 438; % # of rows
% clmz = 566; % # of columns. Shape is rwz x clmz
% lon_flat     = reshape(lon, [rwz*clmz,1]); % flatten lon array
% lat_flat     = reshape(lat, [rwz*clmz,1]); % flatten lat array
% lonlat_flat  = [lon_flat, lat_flat];      % bring them together
% coord_target = [-48.82830, 67.06650];     % kanm aws location
% 
% % see which lat/lon pair is closest to my given coordinates:
% euclidean_from_target = sqrt((lonlat_flat(:,1) - coord_target(1)).^2 +...
%     (lonlat_flat(:,2) - coord_target(2)).^2);
% 
% % id the row in which we have those closest coordinates:
% indx_min_euc = find(euclidean_from_target == min(euclidean_from_target));
% 
% % see what are the coord and make sure that they're close enough:
% targ_coord = lonlat_flat(indx_min_euc,:); 
% 
% indx_col_target = ceil(indx_min_euc / rwz); % find col that has these coords
% indx_row_target = indx_min_euc - (rwz * (indx_col_target - 1)); % find row
% addr_in_template = [indx_row_target, indx_col_target]; % row col in LAT LON
% racmo_swsd = swsd(indx_row_target, indx_col_target, :, :); % get alb from our coords!
% rracmo(366,1) = NaN; % create empty array
% rracmo(:,1)   = racmo_swsd(:,:,:,:); % fill it up!
% 
% racmo_swsd_daily = rracmo; % rename for better comprehension
% 
% 
% 
% 
% %% 3) turn daily to hourly
% 
% % a) turn it from daily to hourly ***** ***** *****
% Nfiles = length(racmo_swsd_daily); % how many days (leap or not)
% hrs    = 24; % hours in a day
% RS((Nfiles*hrs),1) = NaN; % prep hourly array
% RS(1:(length(RS)),1) = NaN; % turn zeros to NaN
% 
% % b) fill every day's 10:00 (am) row ***** ***** *****
% for i = 1:Nfiles % for each day
% % NOTE: time will start w/ Jan1 00:00, so to put for 10:00, use row 11
%     mar_row = 11; % data collected at 10:30 so will put as 10am
%     indx      = mar_row + (hrs*(i-1)); % which row
%     RS(indx,1) = racmo_swsd_daily(i,1); % insert value
% end
% 
% 
% 
% % c) extrapolate ***** ***** *****
% z        = RS;
% z0       = z';  % flip rows and columns for horzcat to work
% z1       = ~isnan(z0);  % see which columns have actual values
% z2       = cumsum(z1-diff([1,z1])/2); % make some magic
% full_alb = interp1(1:nnz(z1),z(z1),z2); % interpolate
% % fill in gaps in front and end 
% %---------front:
% G1 = find(~isnan(full_alb),1); % determine col that has 1st non-Nan value
% G2 = full_alb(1,G1); % extract the actual value of the first non-Nan datapoint
% full_alb(1,1:G1-1) = G2; % fill the NaN gap with the 1st non-Nan value
% %---------back:
% E1 = find(~isnan(full_alb),1,'last'); % position that has last non-Nan
% E2 = full_alb(1,E1); % get very last available non-Nan value
% full_alb(1,E1+1:end) = E2; % fill in last gap in the end
% % revenir à l'original
% RS = full_alb'; % update mar_alb
% 
% clear z z0 z1 z2 full_alb G1 G2 E1 E2 indx i Nfiles mar_row hrs;
% 
% 
% 
% 
% %% 4) time 
% % daily (366)
% load('/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/drive/met.mat');
% timez = met.Time; clear met;
% ceres.Time = timez; ceres.swd = cerez;
% ceres = struct2table(ceres);
% ceres_swd_hourly_adj = table2timetable(ceres);
% ceres_swd_daily_adj_2016 = cerez;
% 
% 
% 
% 
% %% 5) plot
% x = (1:366)';
% y1 = ceres_swd_daily_adj_2016;
% y2 = racmo_swsd_daily;
% y3 = mar_swd_daily;
% % plot(x,y1,'LineWidth',2.5,'DisplayName','CERES','Color',"#77AC30")
% hold on
% % plot(x,y2,'LineWidth',2.5,'DisplayName','RACMO','Color',"#A2142F")
% hold on;
% % plot(x,y1,'LineWidth',2.5,'DisplayName','CERES','Color',"#77AC30")
% 
% plot(x,y3,'LineWidth',2.5,'DisplayName','MAR','Color',"#EDB120")
% hold on
% plot(x,y2,'LineWidth',2.5,'DisplayName','RACMO','Color',"#A2142F")
% 
% 
% 
% set(gcf,'Position',[200, 200, 1500, 1000]) % controls figure size
% fntsz = 15; % font size
% ax = gca;  %get variable to control axes
% fnt = 'Helvetica'; % font type; e.g., times, arial, helvetica, damascus
% title('Comparison of daily SWD for MAR and RACMO for 2016')
% ylabel('Shortwave Down (W/m-2)','fontsize',fntsz,'fontname',fnt); 
% xlabel('Time','FontSize',fntsz,'fontname',fnt);
%     set(gca,'fontsize',fntsz); % enlargen size of axis numbers
%     ax.TitleFontSizeMultiplier = 1.5; % make title larger
%     ax.LineWidth = 1.5; grid on;
% lgd = legend; lgd.FontSize = 16; lgd.Location = 'east';