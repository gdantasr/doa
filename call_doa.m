% Reads data from INMETRO audio recordings and runs a DOA estimation for
% each test. Don't forget to check simulation parameters before running:
%
% micsID:   Select 2 microphones from the circular array (1-11)
% fs:       Sampling frequency for DOA estimation
% N:        Block size for DOA estimation
% carID:    Pass by car(s)
% speed:    Pass by speed(s)

%clear all; 
clc; 
close all;

%% Define parameters and load files

% Data paths
addpath (genpath('../data/'))

% Parameters definition
load('array_circular_11mics.mat');          % Array geometry
load('wb_real');                            % Wheelbase distances
D = load('dist.txt');                       % Car-array distances
array = posMic;
vs = 343;                                   % Sound speed propagation
originalFs = 25600;                         % From acquisition
fs = 25600;                                 % For processing
plot_doa = false;                           % DOA plot flag
plot_fit = false;                           % Fitted curves plot flag
pre_filter = false;                         % Flag to filter signals before processing
micsID = [4 10] ;                            % Mics. for DOA algorithms
methodNames = {'gcc', 'itd', 'lms', 'evd', 'mvdr'}; 

% Calculate distance between two mics. (for DOA)
p1 = array (micsID(1), :);
p2 = array (micsID(2), :);
d = norm (p1-p2);  % Euclidean distance between two points
dist_mic_mic = d;
% Distances structure
%dist.source_mic = 2;
dist.mic_mic = d;
dist.heigth = 1.2;
%dist.wheelbase = 2.5; % Average front-rear axles distance (m)

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
   n90{fileID} = calc_t90(normData(:,1), fs);                           % Azimuth = 90� instant calculation
   cropData{fileID} = originalData{fileID}(n90{fileID}-5*fs:n90{fileID}+5*fs, :);   % Crops data around n90
   
   % Get speed
   v{fileID} = eval([fileNames{fileID}, '.speed;']);                   % Vehicle speed
   
   % Get orientation
   fromto{fileID} = eval([fileNames{fileID}, '.orientation;']);        % Passby orientation
end

%% DOA Estimation
mean_error = zeros(length(fileNames) , length(methodNames)); % Store error for each test
wheelbase = zeros(length(fileNames), 1);
wb_estimated = zeros(length(fileNames), 5);
v_estimated = zeros(length(fileNames), 5);
for fileID = 1 : length(fileNames) % For each file
    
    % Pre processing
    if pre_filter
        fm = 200;       % Low cutoff freq
        fc = 3000;      % High cutoff freq
        [b,a] = firpm(100, [0 fm/(originalFs/2) fc/(originalFs/2) 1], [0 1 1 0]);   % Bandpass filter
%         [b,a] = butter(6, [fm/(originalFs/2) fc/(originalFs/2)], 'bandpass');   % Bandpass filter
    %         [b,a] = butter(5, fc/(originalFs/2), 'low');   % Lowpass filter
    %         [b,a] = butter(5, fm/(originalFs/2), 'high');   % Highpass filter 
        data = filter(b, a, cropData{fileID});
        fmax_itd = min([fs/2, fc]);
    else
        data = cropData{fileID};
        fmax_itd = fs/2;
    end
    
    % Find window size
    N = 2^(ceil(log2( max_samp_triang(5, v{fileID}, fs) )));
        
    % DoA Algorithms
    for methodID = 1:4

        switch methodNames{methodID}           
            case 'gcc'
                [phi, tdd, t, tau, Cmat] = doa_gcc(data(:, micsID(1)), data(:, micsID(2)), d, N, fs);
                th = 0.1;
                
            case 'itd'
                [phi, tdd, t, tau, Cmat] = doa_itd(data(:, micsID(1)), data(:, micsID(2)), d, 4*N, fs, fmax_itd);
                th = [];
                
            case 'lms'
                [phi, tdd, t, tau, Cmat] = doa_lms(data(:, micsID(1)), data(:, micsID(2)), d, N/4, N, 0.25, fs);
                th = 0.1;
            
            case 'evd'
                [phi, tdd, t, tau, Cmat] = doa_aevd(data(:, micsID(1)), data(:, micsID(2)), d, N/4, N, 0.25, fs);
                th = 0.1;
                
            case 'mvdr'
                normData = [data(:, micsID(1))/max(data(:, micsID(1))) data(:, micsID(2))/max(data(:, micsID(2)))];
                [phi, t, tau, Cmat] = doa_mvdr(normData(:, 1), normData(:, 2), d, N, fs);
                th = 0.4; 
                
        end

        % Azimuth = 90� instant calculation
        t90 = t(round(length(t)/2)); % Central point 
        
        % Get speed curve
        vCurve = track_speed(fileNames{fileID}, t, t90);
        is_zero = find(~vCurve);
        vCurvePlot = vCurve;            % Keep zeros to plot
        vCurve(is_zero) = v{fileID};    % Remove zeros

        % Theoretical DOA
        
        % Verifica dist�ncia que o carro passou
        dist.source_mic = D(fileID);
        
        % Verifica sentido de movimento do ve�culo
        if strcmp (fromto{fileID}, '0� -> 180�')
            sinal = -1; 
        else
            sinal = 1;
        end
        carID = strfind('bfjm', fileNames{fileID}(1));  % Check the vehicle from current test
        dist.wheelbase = wb_real(carID);                % Check wheelbase from current vehicle
        t90_inf = t90 + dist.wheelbase / (2*v{fileID}/3.6);
        t90_sup = t90 - dist.wheelbase / (2*v{fileID}/3.6);
        tdd_inf_theo = passby(t90_inf, dist.source_mic, vCurve, t, dist_mic_mic, dist.heigth, sinal);
        tdd_sup_theo = passby(t90_sup, dist.source_mic, vCurve, t, dist_mic_mic, dist.heigth, sinal);
        tdd_inf_theo2 = passby(t90_inf, 2, vCurve, t, dist_mic_mic, dist.heigth, sinal);
        tdd_sup_theo2 = passby(t90_sup, 2, vCurve, t, dist_mic_mic, dist.heigth, sinal);
        
        % Post processing
        [tdd_inf, tdd_sup, fit_inf, fit_sup] = fit_curve2mat (t, t90, tau, Cmat, v{fileID}, fromto{fileID}, plot_fit, 'exclude', [methodNames{methodID},'_',fileNames{fileID}], dist, th);       
        phi_inf = 180/pi*real(acos(vs/d*tdd_inf));
        phi_sup = 180/pi*real(acos(vs/d*tdd_sup));
        

    %% Results

        % Error calculation
        fs_tdd = round(1/mean(diff(t)));        % TDD estimation sample rate
        [~, n90_tdd] = min(abs(t - t90));       % n90 sampled to TDD rate
        time_window = (n90_tdd - fs_tdd : n90_tdd + fs_tdd);
        mean_error_sup = mean(abs(tdd_sup(time_window) - tdd_sup_theo(time_window).'));
        mean_error_inf = mean(abs(tdd_inf(time_window) - tdd_inf_theo(time_window).'));
        mean_error(fileID, methodID) = (mean_error_sup + mean_error_inf)/2;
        
        [~, iMinInf] = min(abs(tdd_inf));
        [~, iMinSup] = min(abs(tdd_sup));
        wb_estimated(fileID, methodID) = abs(t(iMinInf) - t(iMinSup))*v{fileID}/3.6;  % Estimated wheelbase distances
%         wheelbase(fileID) = abs(fit_inf.a - fit_sup.a)*v{fileID}/3.6;  % Estimated wheelbase distances
        v_estimated(fileID, methodID) = fit_sup.v;

        % Plot DOA results
        carID = strfind('bfjm', fileNames{fileID}(1));
        titulo =  [upper(methodNames{methodID}),' Algorithm. Speed ', num2str(v{fileID}),' km/h. Car ', num2str(carID),'.'];
        plotRange = (round(2.5*fs_tdd) : round(7.5*fs_tdd));
        tPlot = t(1:length(plotRange));
        plot_tdd_theo(tPlot, tau, [tdd_sup(plotRange) tdd_inf(plotRange), tdd_sup_theo(plotRange)' tdd_inf_theo(plotRange)'], Cmat(:, plotRange), false, fileNames{fileID}, titulo);
%         set(gcf,'Visible', 'on');

        % Saving
        F = getframe(gcf);
        filePath = ['../Pesquisa/Matlab/Resultados/0821/'];
        figName = [filePath, num2str(fileID), '_doa_', methodNames{methodID}, '_fullband_d', round(num2str(1000*d))];
        imwrite(F.cdata, [figName, '.png'], 'png');     % Save .png
%         savefig([figName, '.fig']);                     % Save .fig
         
    end
    
%     wb_estimated(fileID, 5) = abs(t(iMinInf) - t(iMinSup))*v{fileID}/3.6;  % Estimated wheelbase distances
      
    fileID
end

% error_tab = table( mean_error(:,1), mean_error(:,2), mean_error(:,3), mean_error(:,4), 'VariableNames', methodNames, 'RowNames', fileNames);
% writetable(error_tab, [filePath, 'mean_error_2.csv']);

%save([filePath,'wb_estimated'], 'wb_estimated')
%save([filePath,'v_estimated'], 'v_estimated')

% error_array = table2array(error_tab);
% 
% close all
% figure
% hold on
% for i = 1:size(error_array,2)
%     plot(error_array(:,i), '-o')
% end
% legend(methodNames)
% xlabel('Teste')
% ylabel('Erro m�dio')
% F = getframe(gcf);
% imwrite(F.cdata, [filePath, 'mean_error_2.png'], 'png');  % Save .png