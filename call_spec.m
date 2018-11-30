%   Spectrogram plots for pass by recordings
%   Full band signals

close all; clear; clc

%% Set parameters and load files

fs = 25600;     % Audio sampling frequency
N = 256;        % Spectrogram window size
                    
% Data path
addpath (genpath('../Dissertação/Matlab/data/'))

% Select pass by tests from avaible data

% Pass by car, speed and index
carID = {'f'};                      % options: {'b', 'f', 'j', 'm'};
speed = {'30','50','70','_ac'};     % options: {'30', '50', '60', '70', '_ac'};
passbyID = {'_1'};                  % options: {'1','2','3'}

% Possible file names for specified car and speed
fileNames = repmat(carID, [length(speed), 1]);
fileNames = fileNames(:);
fileNames = strcat( fileNames, repmat(speed', [length(carID), 1] ) );
fileNames = repmat(fileNames', [length(passbyID)],1);
fileNames = fileNames(:);
fileNames = strcat( fileNames, repmat(passbyID', [length(carID)*length(speed), 1]) );

% Check if possible files exist
allFileNames = what(('../Dissertação/Matlab/data/estruturas')); % All files found in dir
allFileNames = allFileNames.mat;    % Only .mat files
for k=1:length(allFileNames)        % Ignores end of strings ('.mat')
    allFileNames{k} = allFileNames{k}(1:end-4);
end
i = ismember(fileNames, allFileNames);
fileNames = fileNames(i);

% Load files
data = cell(length(fileNames),1);
cropData = cell(length(fileNames),1);
for fileID = 1 : length(fileNames)
   load(fileNames{fileID});     
   data{fileID} = eval([fileNames{fileID}, '.soundPressure']);  % Stores audio data
   data{fileID} = data{fileID} - repmat(mean(data{fileID}), [length(data{fileID}), 1] );    % Corrects offset
   n90 = calc_t90(data{fileID}(:,1), fs);                       % Azimuth = 90° instant calculation
   cropData{fileID} = data{fileID}(n90-3*fs:n90+3*fs, 1);       % Crops data around n90
   %data{fileID} = data{fileID} / max(data{fileID});            % Normalization
end

%% Plot spectrogram

% Plot options
plotOption = 3;     % 0 - same test, all sensors (11)
                    % 1 - same sensor, single const. speed vs acceleration fulltime
                    % 2 - same sensor, single const. speed vs acceleration cropped
                    % 3 - same sensor, all const. speed vs acceleration fulltime
                    % 4 - same sensor, all const. speed vs acceleration cropped

switch plotOption
    case 0
        figName = ['../Dissertação/Matlab/plots/spec_allsensors_', fileNames{fileID}];
    case 1 
        figName = ['../Dissertação/Matlab/plots/spec_', fileNames{1},'_fulltime'];
    case 2 
        figName = ['../Dissertação/Matlab/plots/spec_', fileNames{1},'_cropped']; 
        data = cropData;
    case 3 
        figName = ['../Dissertação/Matlab/plots/spec_', carID{1}, '_fulltime'];
    case 4 
        figName = ['../Dissertação/Matlab/plots/spec_', carID{1}, '_cropped']; 
        data = cropData;
end

F = plot_spec(data, fs, N, speed, carID, plotOption);
imwrite(F.cdata, [figName,'.png'], 'png');  % Save .png
savefig([figName, '.fig']);                 % Save .fig