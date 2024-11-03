% rename nc files
% purpose: change names in nc files of merra2 rad for albedo; so that all can be in order
% created on july 15, 2023; Maxim Altan-Lu Shapovalov
% modified for M2T1NXSLV (2m air temp) on oct 29, 2023


% Table of Contents (sections):
% I) making it happen (rename nc files and put in separate folder for check)
% II) seeing which file is extra (actually shows what files you're missing)


% Order of operation:
% 1) Before: went to GES DISC (NASA), downloaded MERRA-2 data using the terminal
    % (read the .txt file), save it in the input folder
% 2) fix ALL the paths here and make sure you're tapping into the right folders
% 3) Run section I and make sure to keep track of the 'met' variable
% 4) Manually check in the created 'ordered' folder that all is good
% 5) Transfer files from 'ordered' to the original folder
% 6) Double check that number adds up to the total number of files in other met files
% 7) Repeat steps 3-6 for other mets
% 8) IF you see that met files have different numbers of files, run section II
    % to see which files exist in one folder and not the other. Hopefully, there
    % aren't that many files missing and you can download them quickly off the GES
    % DISC site again, read in terminal, rename manually if necessary, and plug into
    % the folder that was missing that day
%) 9) After: use the "read_MERRA2_glc.m" to actually get the data



%% 
clear
clc
mp = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/'; % mp = max's path


% I) making it happen

% which met
% met = 'KANL';
% met = 'KANM';
% met = 'KANU';

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

met = aws; % rename to keep it consisent for the existing code

% 1) 401s ----------------------------------------------------------------------
% Define the folder path containing the .nc files

% folderPath = [mp 'drive_python/input/MERRA2/M2T3NXGLC_5.12.4/tair/' met '/'];
folderPath = ['/Users/maxims/Downloads/2. copy/' met '/'];

% Get a list of .nc files with "401" in the filename
filePattern = fullfile(folderPath, 'MERRA2_401*.nc');
fileList = dir(filePattern);

% Create a new folder to save the renamed files
% newFolderPath = [mp 'drive_python/input/MERRA2/M2T3NXGLC_5.12.4/ordered/' met '/'];
newFolderPath = ['/Users/maxims/Downloads/3. renamed/' met '/'];

% mkdir(newFolderPath); % create folder in case it hasn't been created already

% Loop through the files and rename them
for i = 1:numel(fileList)
    % Get the current file name
    oldFileName = fullfile(folderPath, fileList(i).name);
    
    % Create the new file name by replacing "401" with "400"
    newFileName = strrep(fileList(i).name, 'MERRA2_401', 'MERRA2_400');
    
    % Construct the full paths for the old and new files
%     oldFilePath = fullfile(folderPath, oldFileName);
    oldFilePath = oldFileName;
    newFilePath = fullfile(newFolderPath, newFileName);
    
    % Rename the file
    movefile(oldFilePath, newFilePath);
end

% 2) 300s ----------------------------------------------------------------------
% Define the folder path containing the .nc files
% folderPath = [mp 'drive_python/input/MERRA2/M2T3NXGLC_5.12.4/tair/' met '/'];
folderPath = ['/Users/maxims/Downloads/2. copy/' met '/'];

% Get a list of .nc files with "300" in the filename
filePattern = fullfile(folderPath, 'MERRA2_300*.nc');
fileList = dir(filePattern);

% Create a new folder to save the renamed files
% newFolderPath = [mp 'drive_python/input/MERRA2/M2T3NXGLC_5.12.4/ordered/' met '/'];
newFolderPath = ['/Users/maxims/Downloads/3. renamed/' met '/'];

% Loop through the files and rename them
for i = 1:numel(fileList)
    % Get the current file name
    oldFileName = fullfile(folderPath, fileList(i).name);
    
    % Create the new file name by replacing "300" with "400"
    newFileName = strrep(fileList(i).name, 'MERRA2_300', 'MERRA2_400');
    
    % Construct the full paths for the old and new files
%     oldFilePath = fullfile(folderPath, oldFileName);
    oldFilePath = oldFileName;
    newFilePath = fullfile(newFolderPath, newFileName);
    
    % Rename the file
    movefile(oldFilePath, newFilePath);
end


%% II) seeing which file is extra
mp = '/Users/maxims/Library/CloudStorage/OneDrive-UniversityOfOregon/Lab/matlab/icemodel/'; % mp = max's path

folder1 = [mp 'drive_python/input/MERRA2/M2T1NXSLV_5.12.4/KANM/'];  % Specify the path to the first folder
folder2 = [mp 'drive_python/input/MERRA2/M2T1NXSLV_5.12.4/KANL/']; % Specify the path to the second folder

files1 = dir(fullfile(folder1, '*.nc'));  % Get a list of files in folder1
files2 = dir(fullfile(folder2, '*.nc'));  % Get a list of files in folder2

numFiles1 = numel(files1);  % Number of files in folder1
numFiles2 = numel(files2);  % Number of files in folder2

if numFiles1 > numFiles2
    % Folder1 has more files
    extraFiles = setdiff({files1.name}, {files2.name});
    disp('Extra file(s) found in folder1:');
    for i = 1:numel(extraFiles)
        fprintf('%s\n', extraFiles{i});
    end
elseif numFiles2 > numFiles1
    % Folder2 has more files
    extraFiles = setdiff({files2.name}, {files1.name});
    disp('Extra file(s) found in folder2:');
    for i = 1:numel(extraFiles)
        fprintf('%s\n', extraFiles{i});
    end
else
    disp('The two folders have the same number of files.');
end

