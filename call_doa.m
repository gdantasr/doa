clear all; clc; close all;


% Parameters definition
load('array_circular_11mics.mat');          % Array geometry
array = posMic;
vs = 340;                                   % Sound speed propagation
micsID = [4 10] ;                           % Mics. for DOA algorithms
originalFs = 25600;                         % From acquisition
fs = 32000;                                 % For processing
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

% Data paths
addpath ('/home/netware/users/gabidantas/Documents/Mestrado/Dissertação/Matlab/data/')
addpath ('estruturas');     % Path to NEW data (INMETRO recordings)


% Store all file names possibilities to go through
carID = {'b', 'f', 'j', 'm'};               % Vehicles ID 
speed = {'30', '50', '60', '70', '_ac'};    % Speeds (estimated)
allFileNames = repmat(carID, [length(speed), 1]); 
allFileNames = allFileNames(:);
allFileNames = strcat( allFileNames, repmat(speed', [length(carID), 1] ) );

for namesID = 1 : length(allFileNames) % For each possible filename...
    
    passbyID = 1;
    filename = allFileNames{namesID};
    filename = strcat( filename, '_', num2str(passbyID) );
    
    % Check if possible name is actual name
    while exist(strcat(filename, '.mat'), 'file') == 2
        
        % Load file and store relevant data
        load(filename)
        eval(['v = ' , filename, '.speed;'])                        % Vehicle speed
        eval(['fromto = ' , filename, '.orientation;'])             % Passby orientation
        eval(['originalData = ', filename, '.soundPressure;'])   % Sound pressure measure
       
        % GAMBIARRA enquanto não conserta o angulo nos structs
        if fromto == '0 -> 180'
            fromto = '180 -> 0';
        elseif fromto == '180 -> 0'
            fromto = '0 -> 180';
        end
        
        %% DOA Estimation
        
        % Pre processing
        fm = 200;       % Low cutoff freq
        fc = 4000;      % High cutoff freq
        [b,a] = butter(5, [fm/(originalFs/2) fc/(originalFs/2)], 'bandpass');   % Bandpass filter 
        filtData = filter(b, a, originalData);     
        resampData = resample(filtData, fs, originalFs);
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
        [tdd_inf, tdd_sup, fit_inf, fit_sup] = edge_detect (t, t90, tau, Cmat, v, fromto, plot_fit, 'exclude', filename, dist);       
        phi_inf = 180/pi*real(acos(vs/d*tdd_inf)); 
        phi_sup = 180/pi*real(acos(vs/d*tdd_sup));
        
        % Plot DOA results
        titulo =  ['GCC-PHAT function. TestID: ', filename];
        plot_tdd(t, tau, [tdd_sup tdd_inf], Cmat, false, filename, titulo);
        pause
        
        % Filename update
        passbyID = passbyID + 1;
        filename = allFileNames{namesID};
        filename = strcat( filename, '_', num2str(passbyID) );
        
    end
    
end