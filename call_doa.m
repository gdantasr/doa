% Reads data from INMETRO audio recordings and runs a DOA estimation for
% each test. Don't forget to check simulation parameters before running:
%
% micsID:   Select 2 microphones from the circular array (1-11)
% fs:       Sampling frequency for DOA estimation
% N:        Block size for DOA estimation
% carID:    Pass by car(s)
% speed:    Pass by speed(s)

% clear all; 
clc; 
close all;

%% Define parameters and load files

% Data paths
addpath (genpath('../data/'))

% Parameters definition
load('array_circular_11mics.mat');          % Array geometry
array = posMic;
vs = 340;                                   % Sound speed propagation
originalFs = 25600;                         % From acquisition
fs = 25600;                                 % For processing
plot_doa = false;                           % DOA plot flag
plot_fit = false;                           % Fitted curves plot flag 
micsID = [5 9] ;                           % Mics. for DOA algorithms

% Calculate distance between two mics. (for DOA)
p1 = array (micsID(1), :);
p2 = array (micsID(2), :);
d = norm (p1-p2);  % Euclidean distance between two points
dist_mic_mic = d;
% Distances structure
dist.source_mic = 2;
dist.mic_mic = d;
dist.heigth = 1.2;
dist.axle = 2.5; % Average front-rear axles distance (m)

% Experiment distances in meters (INMETRO)
dist_source_mic = 2;
heigth_mic = 1;

% Possible file names for specified car and speed
carID = {'b', 'f', 'j', 'm'};               % Vehicles ID 
speed = {'30', '50', '60', '70', '_ac'};    % Speeds (estimated)
passbyID = {'_1','_2','_3'};                   
fileNames = repmat(carID, [length(speed), 1]);
fileNames = fileNames(:);
fileNames = strcat( fileNames, repmat(speed', [length(carID), 1] ) );
fileNames = repmat(fileNames', [length(passbyID)],1);
fileNames = fileNames(:);
fileNames = strcat( fileNames, repmat(passbyID', [length(carID)*length(speed), 1]) );

% Check if possible files exist
allFileNames = what(('../data/estruturas'));    % All files found in dir
allFileNames = allFileNames.mat;                % Only .mat files
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
n90 = cell(length(fileNames),1);
for fileID = 1 : length(fileNames)
   load(fileNames{fileID});     % Load struct
   
   % Get audio data
   originalData{fileID} = eval([fileNames{fileID}, '.soundPressure']);  % Stores audio data
   bias = repmat(mean(originalData{fileID}), [length(originalData{fileID}), 1]);
   originalData{fileID} = originalData{fileID} - bias;                  % Remove bias
   normData = originalData{fileID} / max(abs(originalData{fileID}(:)));
   n90{fileID} = calc_t90(normData(:,1), fs);                           % Azimuth = 90° instant calculation
   cropData{fileID} = originalData{fileID}(n90{fileID}-3*fs:n90{fileID}+3*fs, 1);   % Crops data around n90
   
   % Get speed
   v{fileID} = eval([fileNames{fileID}, '.speed;']);                   % Vehicle speed
   
   % Get orientation
   fromto{fileID} = eval([fileNames{fileID}, '.orientation;']);        % Passby orientation
end


%% DOA Estimation
mean_error = zeros(length(fileNames) , 1); % Store error for each test
for fileID = 1 : length(fileNames) % For each file
    
    % Pre processing
    fm = 200;       % Low cutoff freq
    fc = 3000;      % High cutoff freq
    [b,a] = butter(5, [fm/(originalFs/2) fc/(originalFs/2)], 'bandpass');   % Bandpass filter
%         [b,a] = butter(5, fc/(originalFs/2), 'low');   % Lowpass filter
%         [b,a] = butter(5, fm/(originalFs/2), 'high');   % Highpass filter 
    filtData = filter(b, a, originalData{fileID});     
    resampData = resample(filtData, fs, originalFs);
%     resampData = resample(originalData{fileID}, fs, originalFs);
    data = resampData;
    N = 2^(ceil(log2( max_samp_triang(5, v{fileID}, fs) )));

    % Azimuth = 90° instant calculation
    t90 = n90{fileID}/fs;

    % DoA Algorithms
%     [phi, tdd, t, tau, Cmat] = doa_gcc(data(:, micsID(1)), data(:, micsID(2)), d, N, fs);
    [phi, tdd, t, tau, Cmat] = doa_itd(data(:, micsID(1)), data(:, micsID(2)), d, 4*N, fs);
%     [phi, tdd, t, tau, Cmat] = doa_lms(data(:, micsID(1)), data(:, micsID(2)), d, N/4, N, 0.25, fs);
%     [phi, tdd, t, tau, Cmat] = doa_aevd(data(:, micsID(1)), data(:, micsID(2)), d, N/4, N, 0.25, fs);

    % Get speed curve
    vCurve = track_speed(fileNames{fileID}, t, t90);
    is_zero = find(~vCurve);
    vCurve(is_zero) = v{fileID}; % Remove zeros
    
    % Theoretical DOA
    % Verifica sentido de movimento do veículo
    if strcmp (fromto{fileID}, '0° -> 180°')
        sinal = -1; 
    else
        sinal = 1;
    end
    t90_inf = t90 + dist.axle / (2*v{fileID}/3.6);
    t90_sup = t90 - dist.axle / (2*v{fileID}/3.6);
    tdd_inf_theo = passby(t90_inf, dist.source_mic, vCurve, t, dist_mic_mic, heigth_mic, sinal);
    tdd_sup_theo = passby(t90_sup, dist.source_mic, vCurve, t, dist_mic_mic, heigth_mic, sinal);       

    % Post processing
    [tdd_inf, tdd_sup, fit_inf, fit_sup] = fit_curve2mat (t, t90, tau, Cmat, v{fileID}, fromto{fileID}, plot_fit, 'exclude', fileNames{fileID}, dist);       
    phi_inf = 180/pi*real(acos(vs/d*tdd_inf));
    phi_sup = 180/pi*real(acos(vs/d*tdd_sup));

%% Results

    % Error calculation
    fs_tdd = round(1/mean(diff(t)));        % TDD estimation sample rate
    [~, n90_tdd] = min(abs(t - t90));       % n90 sampled to TDD rate
    time_window = (n90_tdd - fs_tdd : n90_tdd + fs_tdd);
    mean_error_sup = mean(abs(tdd_sup(time_window) - tdd_sup_theo(time_window).'));
    mean_error_inf = mean(abs(tdd_inf(time_window) - tdd_inf_theo(time_window).'));
    mean_error(fileID) = (mean_error_sup + mean_error_inf)/2;

    % Plot DOA results
    titulo =  ['ITD Algortihm. Fullband signal. Test ID: ', fileNames{fileID}];

%         plot_tdd_speed(t, tau, [tdd_sup tdd_inf], [tdd_sup_theo tdd_inf_theo], Cmat, false, fileNames{fileID}, titulo, vCurve);
    plot_tdd_theo(t, tau, [tdd_sup tdd_inf, tdd_sup_theo' tdd_inf_theo'], Cmat, false, fileNames{fileID}, titulo);
    %set(gcf,'Visible', 'off');

    % Saving
    F = getframe(gcf);
    figName = ['../Dissertação/Matlab/plots/0404/doa_itd_band_', num2str(fm), '_', num2str(fc), '_', fileNames{fileID}, '_d', num2str(100*d)];
    imwrite(F.cdata, [figName,'.png'], 'png');  % Save .png
%         savefig([figName, '.fig']);                 % Save .fig
end
itd_mean_error = mean_error;
