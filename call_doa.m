% Reads data from INMETRO audio recordings and runs a DOA estimation for
% each test. Don't forget to check simulation parameters before running:
%
% micsID:   Select 2 microphones from the circular array (1-11)
% fs:       Sampling frequency for DOA estimation
% N:        Block size for DOA estimation
% carID:    Pass by car(s)
% speed:    Pass by speed(s)

clear all; clc; close all;

%% Define parameters and load files

% Data paths
addpath (genpath('../Dissertação/Matlab/data/'))

% Parameters definition
load('array_circular_11mics.mat');          % Array geometry
array = posMic;
vs = 340;                                   % Sound speed propagation
micsID = [4 10] ;                           % Mics. for DOA algorithms
originalFs = 25600;                         % From acquisition
fs = 25600;                                 % For processing
N = 256;                                    % Block size (samples)
plot_doa = false;                           % DOA plot flag
plot_fit = false;                           % Fitted curves plot flag 

% Calculate distance between two mics. (for DOA)
p1 = array (micsID(1), :);
p2 = array (micsID(2), :);
d = norm (p1-p2);  % Euclidean distance between two points

% Experiment distances in meters (INMETRO)
% TODO : BOTAR NUMA TABELA/STRUCT 
% TODO : GRAVAR FS NO STRUCT
dist_source_mic = 2;
heigth_mic = 1;
dist_mic_mic = d;
% dist_source_source = 2.5; NÃO USA EM LUGAR NENHUM

% Possible file names for specified car and speed
carID = {'j'};      %{'b', 'f', 'j', 'm'};               % Vehicles ID 
speed = {'30'};     %'30', '50', '60', '70', '_ac'};    % Speeds (estimated)
passbyID = {'_1'};   %{'_1','_2','_3'};                   
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
v = cell(length(fileNames),1);
fromto = cell(length(fileNames),1);
for fileID = 1 : length(fileNames)
   load(fileNames{fileID});     % Load struct
   
   % Get audio data
   originalData{fileID} = eval([fileNames{fileID}, '.soundPressure']);  % Stores audio data
   bias = repmat(mean(originalData{fileID}), [length(originalData{fileID}), 1]);
   originalData{fileID} = originalData{fileID} - bias;                  % Remove bias
   n90 = calc_t90(originalData{fileID}(:,1), fs);                       % Azimuth = 90° instant calculation
   cropData{fileID} = originalData{fileID}(n90-3*fs:n90+3*fs, 1);       % Crops data around n90
   
   % Get speed
   v{fileID} = eval([fileNames{fileID}, '.speed;']);                   % Vehicle speed
   
   % Get orientation
   fromto{fileID} = eval([fileNames{fileID}, '.orientation;']);        % Passby orientation
end


%% DOA Estimation
for namesID = 1 : length(fileNames) % For each file
        
    % Pre processing
    fm = 200;       % Low cutoff freq
    fc = 4000;      % High cutoff freq
    [b,a] = butter(5, [fm/(originalFs/2) fc/(originalFs/2)], 'bandpass');   % Bandpass filter
    %[b,a] = butter(5, fc/(originalFs/2), 'low');   % Lowpass filter 
    filtData = filter(b, a, originalData{namesID});     
    resampData = resample(filtData, fs, originalFs);
    %resampData = resample(originalData{namesID}, fs, originalFs);
    data = resampData;

    % Azimuth = 90° instant calculation
    n90 = energy_peak(data(:,1), fs);
    t90 = n90/fs;

    % GCC-PHAT       
    [phi, tdd, t, tau, Cmat] = doa_gcc_mine(data(:, micsID(1)), data(:, micsID(2)), d, N, fs);

    % Post processing
    % TODO - passar distancias como arg (cell? struct?)
    dist.source_mic = 2;
    dist.mic_mic = d;
    dist.heigth = 1;
    [tdd_inf, tdd_sup, fit_inf, fit_sup] = edge_detect (t, t90, tau, Cmat, v{namesID}, fromto{namesID}, plot_fit, 'exclude', fileNames{namesID}, dist);       
    phi_inf = 180/pi*real(acos(vs/d*tdd_inf)); 
    phi_sup = 180/pi*real(acos(vs/d*tdd_sup));

    %fit_sup

    % Plot DOA results
    titulo =  ['GCC-PHAT function. Banda ', num2str(fm), '-', num2str(fc), ' Hz. TestID: ', fileNames{namesID}];
    plot_tdd(t, tau, [tdd_sup tdd_inf], Cmat, false, fileNames{namesID}, titulo);
    F = getframe(gcf);
    figName = ['../Dissertação/Matlab/plots/doa_band_', num2str(fm),'_' , num2str(fc), '_',fileNames{namesID}, '_fs', num2str(fs)];
    imwrite(F.cdata, [figName,'.png'], 'png');  % Save .png
    savefig([figName, '.fig']);                 % Save .fig

end