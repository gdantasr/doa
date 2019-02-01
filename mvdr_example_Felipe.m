clear all
close all
clc



% Create vectors of x- and y-coordinates of microphone positions 

xPos = -0.2:0.2:0.2; % 1xP vector of x-positions [m]
yPos = zeros(1, numel(xPos)); % 1xP vector of y-positions [m]

%yPos = -0.3:0.3:0.3;
%xPos = zeros(1, numel(yPos));

zPos = zeros(1, numel(xPos));
elementWeights = ones(1, numel(xPos))/numel(xPos); % 1xP vector of weightings
%size(elementWeights)

% Define arriving angles and frequency of input signals
%thetaArrivalAngles = [-30 10]; % degrees
%phiArrivalAngles = [0 0]; % degrees
%f = 800; % [Hz]
c = 340; % [m/s]
fs = 44.1e3; % [Hz]

% Define array scanning angles (1D, so phi = 0)
%thetaScanAngles = -90:0.1:90; % degrees
thetaScanAngles = -90:5:90; % degrees
%thetaScanAngles = 180:-5:0; % degrees


phiScanAngles = 0; % degrees


% Create input signal
%inputSignal = createSignal(xPos, yPos, zPos, f, c, fs, thetaArrivalAngles, phiArrivalAngles);
x1 = audioread('E:\COE873\Audios Cortados\teste2_carro1_v30_mic1.wav');
x2 = audioread('E:\COE873\Audios Cortados\teste2_carro1_v30_mic2.wav');
x3 = audioread('E:\COE873\Audios Cortados\teste2_carro1_v30_mic3.wav');




%inputSignal = transpose([x1 x2 x3]);
Xin = [x1 x2 x3];


[Nx,N_D] = size(Xin);

%Nfft = 512;                          % FFT length = time window length
Nfft = 2048;
Nfh = Nfft/2+1;
%M_D = Nfft/4;                          % frame hop size 
M_D = Nfft/2;                          % frame hop size 


%fmax = 4000/(2*0.2);
fmax = 2000;
fmin = 300;
Nfmax = round(Nfft*fmax/fs);
Nfmin = round(Nfft*fmin/fs);

h = hanning(Nfft);
H = h(:) * ones(1,N_D);                % matrix of time windows

%Nfmax

n = 0;

for  m_D = 1:M_D:Nx-Nfft+1

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

[aux,idx]=max(S);
doa = thetaScanAngles(idx);
doa_vec(n) = doa;

spectrumNormalized = abs(S)/max(abs(S));

spectrumLog = 10*log10(spectrumNormalized);
%S_vec(:,n) = spectrumLog(:);
S_vec(:,n) = spectrumNormalized(:);

end

%doa_vec = tsmovavg(doa_vec,'s',10);
doa_vec = medfilt1(doa_vec,5);

X = -90:5:90;
Y = 1:size(S_vec,2);
figure, mesh(Y,X,abs(S_vec));

%size(R)

figure, plot(doa_vec,'o')

pause

%Normalise spectrum
spectrumNormalized = abs(S)/max(abs(S));

%Convert to decibel
spectrumLog = 10*log10(spectrumNormalized);


%Plot array
fig1 = figure;
fig1.Color = 'w';
ax = axes('Parent', fig1);
scatter(ax, xPos, yPos, 20, 'filled')
axis(ax, 'square')
ax.XLim = [-1 1];
ax.YLim = [-1 1];
grid(ax, 'on')
title(ax, 'Microphone positions')

%Plot steered response with indicator lines
fig2 = figure;
fig2.Color = 'w';
ax = axes('Parent', fig2);
plot(ax, thetaScanAngles, spectrumLog)
grid(ax, 'on')
ax.XLim = [thetaScanAngles(1) thetaScanAngles(end)];

% for j=1:numel(thetaArrivalAngles)
%     indx = find(thetaScanAngles >= thetaArrivalAngles(j), 1);
%     line(ax, [thetaScanAngles(indx) thetaScanAngles(indx)], ax.YLim, ...
%         'LineWidth', 1, 'Color', 'r', 'LineStyle', '--');
% end
% xlabel(ax, '\theta')
% ylabel(ax, 'dB')
