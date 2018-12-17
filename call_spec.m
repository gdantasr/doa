%   Spectrogram plots for pass by recordings
%   Full band signals

%close all; clear; clc

%% Set parameters and load files

N = 256;            % Spectrogram window size
originalFs = 25600; % Audio sampling frequency

% Pre processing
fm = 200;       % Low cutoff freq
fc = 4000;      % High cutoff freq
fs = 2*fc;     
%[b,a] = butter(5, [fm/(originalFs/2) fc/(originalFs/2)], 'bandpass');   % Bandpass filter
[b,a] = butter(5, fc/(originalFs/2), 'low');   % Lowpass filter
                    
% Data path
addpath (genpath('../Dissertação/Matlab/data/'))

% Select pass by tests from avaible data

% Pass by car, speed and index
carID = {'j'};                      % options: {'b', 'f', 'j', 'm'};
speed = {'30', '50', '60', '70', '_ac'};     % options: {'30', '50', '60', '70', '_ac'};
passbyID = {'_2'};                  % options: {'1','2','3'}

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
originalData = cell(length(fileNames),1);
cropData = cell(length(fileNames),1);
data = cell(length(fileNames),1);
cropData = cell(length(fileNames),1);
n90 = cell(length(fileNames),1);
filtData = cell(length(fileNames),1);
for fileID = 1 : length(fileNames)
   load(fileNames{fileID});     % Load struct
   
   % Get audio data
   originalData{fileID} = eval([fileNames{fileID}, '.soundPressure']);  % Stores audio data
   bias = repmat(mean(originalData{fileID}), [length(originalData{fileID}), 1]);
   originalData{fileID} = originalData{fileID} - bias;                  % Remove bias
   normData = originalData{fileID} / max(abs(originalData{fileID}(:)));
   n90{fileID} = calc_t90(normData(:,1), originalFs);                           % Azimuth = 90° instant calculation
   cropData{fileID} = originalData{fileID}(n90{fileID}-3*originalFs:n90{fileID}+3*originalFs, 1);   % Crops data around n90
   filtData{fileID} = filter(b, a, originalData{fileID});
   resampData{fileID} = resample(filtData{fileID}, fs, originalFs);
end

%% Plot spectrogram

% Plot options
plotOption = 6;     % 0 - same test, all sensors (11)
                    % 1 - same sensor, single const. speed vs acceleration fulltime
                    % 2 - same sensor, single const. speed vs acceleration cropped
                    % 3 - same sensor, all const. speed vs acceleration fulltime
                    % 4 - same sensor, all const. speed vs acceleration cropped
                    % 5 - same sensor, all const. speed vs acceleration cropped and filtered

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
    case 5
        figName = ['../Dissertação/Matlab/plots/spec_', carID{1}, '_filt_fulltime'];
        data = resampData;
    case 6
        figName = ['../Dissertação/Matlab/plots/spec_', carID{1}, '_mean'];
        data = resampData;    
end

F = plot_spec(data, fs, N, speed, carID, plotOption, n90);
% imwrite(F.cdata, [figName,'.png'], 'png');  % Save .png
% savefig([figName, '.fig']);                 % Save .fig