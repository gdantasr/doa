function [doa_vec, t, tau, S_vec] = doa_mvdr(Xin1, Xin2, d, Nfft, fs)

% clear all
% close all
% clc

Xin(:,1) = Xin1;
Xin(:,2) = Xin2;


yPos = [-d/2 d/2];
xPos = zeros(1, numel(yPos));
zPos = zeros(1, numel(xPos));

elementWeights = ones(1, numel(xPos))/numel(xPos); % 1xP vector of weightings

% Define arriving angles and frequency of input signals
c = 343; % [m/s]

thetaScanAngles = 90; % degrees
%thetaScanAngles = atan(2/1.2)*180/pi;
%phiScanAngles = -90:1:90; % degrees

vs = 343;  
maxDelay = 2 + ceil( d/vs*fs ); 
L=4;
tau = linspace(-maxDelay/fs, maxDelay/fs, L*2*maxDelay);
X=tau;
phiScanAngles = (90-180/pi*real(acos(vs/d*tau)));

[Nx,N_D] = size(Xin);

%Nfft = 512;                          % FFT length = time window length
Nfft = 1024;
Nfh = Nfft/2+1;
%M_D = Nfft/4;                          % frame hop size 
M_D = Nfft/2;                          % frame hop size 

%fmax = 4000/(2*0.2);
fmax = 2000;
fmin = 100;
Nfmax = round(Nfft*fmax/fs);
Nfmin = round(Nfft*fmin/fs);

h = hanning(Nfft);
H = h(:) * ones(1,N_D);                % matrix of time windows

%Nfmax

n = 0;

for  m_D = 1:M_D:Nx-Nfft+1
    
    n = n + 1;

    m1 = m_D:m_D+Nfft-1;

    X_D = fft(Xin(m1,:) .* H).';        % apply FFT to weighted frames of all mics

    inputSignal = transpose(X_D);
    %size(inputSignal)
    % Create steering vector/matrix
    %e = steeringVector(xPos, yPos, zPos, f, c, thetaScanAngles, phiScanAngles);

    S = 0;

    for Nf=Nfmin:Nfmax

        f=fs/Nfft*Nf;
        e = steeringVector(xPos, yPos, zPos, f, c, thetaScanAngles, phiScanAngles);

        % Create cross spectral matrix
        R = crossSpectralMatrix(transpose(inputSignal(Nf,:)));

        % Calculate delay-and-sum steered response
        S = S + steeredResponseMVDR(R, e, elementWeights);

    end

    [aux,idx]=max(S);
    %doa = thetaScanAngles(idx);
    doa = phiScanAngles(idx);
    doa_vec(n) = doa;
    spectrumNormalized = abs(S)/max(abs(S));
    spectrumLog = 10*log10(spectrumNormalized);
    %S_vec(:,n) = spectrumLog(:);
    S_vec(:,n) = spectrumNormalized(:);

end

doa_vec = medfilt1(doa_vec,5);

%X = -90:1:90;
Y = 1:size(S_vec,2);
doa_vec = doa_vec*100/35;

t = M_D/fs*[0:length(Y)-1];