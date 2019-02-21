%clear all; close all; clc;
close all

%% Define parameters and load files

% Data paths
addpath (genpath('../data/'))

% Parameters definition
load('array_circular_11mics.mat');          % Array geometry
array = [posMic zeros(length(posMic), 1)];  % Add z coordinates
%array = posMic;
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


%%
for namesID = 1 : length(fileNames) % For each file

    % Create vectors of weights
    elementWeights = ones(1, length(array)) / length(array); % 1xP vector of weightings

    % Define arriving angles and frequency of input signals
    %thetaArrivalAngles = [-30 10]; % degrees
    %phiArrivalAngles = [0 0]; % degrees
    %f = 800; % [Hz]
    c = 340; % [m/s]
    %fs = 44.1e3; % [Hz]

    % Define array scanning angles (1D, so phi = 0)
    % TODO. MUDAR AQUI
    %thetaScanAngles = -90:0.1:90; % degrees
    %thetaScanAngles = -90:5:90; % degrees
    %thetaScanAngles = 180:-5:0; % degrees 
    %phiScanAngles = zeros(size(thetaScanAngles)); % degrees
    %phiScanAngles = 180;

    phiScanAngles = -90:5:90; % degrees
    thetaScanAngles = 90;
    
    % Input signal
    Xin = originalData{namesID};


    [Nx, N_D] = size(Xin);

    %Nfft = 512;                          % FFT length = time window length
    Nfft = 2048;
    Nfh = Nfft/2+1;
    %M_D = Nfft/4;                          % frame hop size 
    M_D = Nfft/2;                          % frame hop size 


    %fmax = 4000/(2*0.2);
    fmax = 4000;
    fmin = 2000;
    Nfmax = round(Nfft*fmax/fs);
    Nfmin = round(Nfft*fmin/fs);

    h = hanning(Nfft);
    H = h(:) * ones(1,N_D);                % matrix of time windows

    %Nfmax

    n = 0;
    for m_D = 1:M_D:Nx-Nfft+1

        %disp(n)    
        n = n + 1;
        m1 = m_D:m_D+Nfft-1;
        X_D = fft(Xin(m1,:) .* H).';        % apply FFT to weighted frames of all mics
        inputSignal = transpose(X_D);
        %size(inputSignal)
        % Create steering vector/matrix
        %e = steeringVector(xPos, yPos, zPos, f, c, thetaScanAngles, phiScanAngles);

        S = 0;

        %for Nf=1:Nfmax
        for Nf=Nfmin:Nfmax
            f=fs/Nfft*Nf;
            
            %e = steeringVector(xPos, yPos, zPos, f, c, thetaScanAngles, phiScanAngles);
            xPos = array(:, 3);
            yPos = array(:, 1);
            zPos = array(:, 2);

            e = steeringVector(xPos, yPos, zPos, f, c, thetaScanAngles, phiScanAngles); 
            
            %Nfmax
            % Create cross spectral matrix
            R = crossSpectralMatrix(transpose(inputSignal(Nf,:)));
            
            %Nf
            %size(R)
            %pause

            %size(R)

            % Calculate delay-and-sum steered response
            S = S + steeredResponseDelayAndSum(R, e, elementWeights);

        end

        [aux,idx] = max(S);
        doa = phiScanAngles(idx);
        doa_vec(n) = doa;

        spectrumNormalized = abs(S)/max(abs(S));
        spectrumLog = 10*log10(spectrumNormalized);
        %S_vec(:,n) = spectrumLog(:);
        S_vec(:,n) = spectrumNormalized(:);

    end

    %doa_vec = tsmovavg(doa_vec,'s',10);
    doa_vec = medfilt1(doa_vec,5);

    %X = -90:5:90;
    X = phiScanAngles;
    Y = 1:size(S_vec,2);
    figure, mesh(Y,X,abs(S_vec));

    %size(R)

    figure, plot(doa_vec,'o')

    %pause

    %Normalise spectrum
    spectrumNormalized = abs(S)/max(abs(S));

    %Convert to decibel
    spectrumLog = 10*log10(spectrumNormalized);


    %Plot array
    fig1 = figure;
    fig1.Color = 'w';
    ax = axes('Parent', fig1);
    scatter(ax, yPos, zPos, 20, 'filled')
    axis(ax, 'square')
    ax.XLim = [-1 1];
    ax.YLim = [-1 1];
    grid(ax, 'on')
    title(ax, 'Microphone positions')

    %Plot steered response with indicator lines
    fig2 = figure;
    fig2.Color = 'w';
    ax = axes('Parent', fig2);
    plot(ax, phiScanAngles, spectrumLog)
    grid(ax, 'on')
    ax.XLim = [phiScanAngles(1) phiScanAngles(end)];

    % for j=1:numel(thetaArrivalAngles)
    %     indx = find(thetaScanAngles >= thetaArrivalAngles(j), 1);
    %     line(ax, [thetaScanAngles(indx) thetaScanAngles(indx)], ax.YLim, ...
    %         'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
    % end
    % xlabel(ax, '\theta')
    % ylabel(ax, 'dB')

end