% Reads data from INMETRO audio recordings and runs a DOA estimation for
% each test. Don't forget to check simulation parameters before running:
%
% micsID:   Select 2 microphones from the circular array (1-11)
% fs:       Sampling frequency for DOA estimation
% N:        Block size for DOA estimation
% carID:    Pass by car(s)
% speed:    Pass by speed(s)

%clear all; clc; 
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

% Experiment distances in meters (INMETRO)
% TODO : BOTAR NUMA TABELA/STRUCT 
% TODO : GRAVAR FS NO STRUCT
dist_source_mic = 2;
heigth_mic = 1;
%dist_mic_mic = d;
% dist_source_source = 2.5; NÃO USA EM LUGAR NENHUM

% Possible file names for specified car and speed
carID = {'b'};%{'b', 'f', 'j', 'm'};               % Vehicles ID 
speed = {'30'};%{'30', '50', '60', '70', '_ac'};    % Speeds (estimated)
passbyID = {'_1'};%{'_1','_2','_3'};                   
fileNames = repmat(carID, [length(speed), 1]);
fileNames = fileNames(:);
fileNames = strcat( fileNames, repmat(speed', [length(carID), 1] ) );
fileNames = repmat(fileNames', [length(passbyID)],1);
fileNames = fileNames(:);
fileNames = strcat( fileNames, repmat(passbyID', [length(carID)*length(speed), 1]) );

% Check if possible files exist
%allFileNames = what(('../Dissertação/Matlab/data/estruturas')); % All files found in dir
allFileNames = what(('../data/estruturas')); % All files found in dir
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
for namesID = 1 : length(fileNames) % For each file
    
    for mic = 4
        
        micsID = [mic 14-mic] ;               % Mics. for DOA algorithms
        % Calculate distance between two mics. (for DOA)
        p1 = array (micsID(1), :);
        p2 = array (micsID(2), :);
        d = norm (p1-p2);  % Euclidean distance between two points
        dist_mic_mic = d;

        % Pre processing
        fm = 500;       % Low cutoff freq
        fc = 2000;      % High cutoff freq
        [b,a] = butter(5, [fm/(originalFs/2) fc/(originalFs/2)], 'bandpass');   % Bandpass filter
        %[b,a] = butter(5, fc/(originalFs/2), 'low');   % Lowpass filter
        %[b,a] = butter(5, fm/(originalFs/2), 'high');   % Highpass filter 
        filtData = filter(b, a, originalData{namesID});     
        %resampData = resample(filtData, fs, originalFs);
        resampData = resample(originalData{namesID}, fs, originalFs);
        data = resampData;
        N = 2^(ceil(log2( max_samp_triang(5, v{namesID}, fs) )));

        % Azimuth = 90° instant calculation
        t90 = n90{namesID}/fs;
      
        % DoA Algorithms
        [phi, tdd, t, tau, Cmat] = doa_gcc(data(:, micsID(1)), data(:, micsID(2)), d, N, fs);
%        [phi, tdd, t, tau, Cmat] = doa_itd(data(:, micsID(1)), data(:, micsID(2)), d, N, fs);
%        [phi, tdd, t, tau, Cmat] = doa_lms(data(:, micsID(1)), data(:, micsID(2)), d, N/4, N, 0.25, fs);
%        [phi, tdd, t, tau, Cmat] = doa_evd(data(:, micsID(1)), data(:, micsID(2)), d, N, fs);
        
        % Post processing
        dist.source_mic = 2;
        dist.mic_mic = d;
        dist.heigth = 1;
        [tdd_inf, tdd_sup, fit_inf, fit_sup] = edge_detect (t, t90, tau, Cmat, v{namesID}, fromto{namesID}, plot_fit, 'exclude', fileNames{namesID}, dist);       
        phi_inf = 180/pi*real(acos(vs/d*tdd_inf)); 
        phi_sup = 180/pi*real(acos(vs/d*tdd_sup));

        % Get speed curve
        vCurve = track_speed(fileNames{namesID}, t, t90);
        
%% Theoretical DOA

        axleDistance = 2.5; % Front-rear axles distance (m)
        % Verifica sentido de movimento do veículo
        if strcmp (fromto{namesID}, '0° -> 180°')
            sinal = -1; 
        else
            sinal = 1;
        end
        t90_inf = t90 + axleDistance/v{namesID};
        t90_sup = t90 - axleDistance/v{namesID};
        tdd_inf_theo = passby( t90_inf, dist.source_mic, vCurve, t, dist_mic_mic, heigth_mic, sinal);
        tdd_sup_theo = passby( t90_sup, dist.source_mic, vCurve, t, dist_mic_mic, heigth_mic, sinal);

%% Results
        
        % Error calculation
        itd_mean_error_sup = mean(abs(tdd_sup-tdd_sup_theo.'));
        itd_mean_error_inf = mean(abs(tdd_inf-tdd_inf_theo.'));
        
        % Plot DOA results
        titulo =  ['Função GCC-PHAT. Banda ', num2str(fm), '-', num2str(fc), ' Hz. TestID: ', fileNames{namesID}];
               
%         plot_tdd_speed(t, tau, [tdd_sup tdd_inf], [tdd_sup_theo tdd_inf_theo], Cmat, false, fileNames{namesID}, titulo, vCurve);
        plot_tdd_speed_theo(t, tau, [tdd_sup tdd_inf, tdd_sup_theo' tdd_inf_theo'], Cmat, false, fileNames{namesID}, titulo, vCurve);
        
        %set(gcf,'Visible', 'off');
        % Saving
        F = getframe(gcf);
        figName = ['../Dissertação/Matlab/plots/doa_method/doa_speed_band_', num2str(fm),'_' , num2str(fc), '_',fileNames{namesID}, '_d', num2str(100*d)];
%         imwrite(F.cdata, [figName,'.png'], 'png');  % Save .png
%         savefig([figName, '.fig']);                 % Save .fig

    end
    
%     % Paste figures on the subplots
%     figure;
%     for i=1:4
%         h(i) = subplot(4,1,i);
%         copyobj(allchild(get(sub(i),'CurrentAxes')),h(i));
%     end

end
