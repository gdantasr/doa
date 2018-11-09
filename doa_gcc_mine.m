function [azimuth, delay, t, tau, Cmat] = doa_gcc(x1, x2, dx, N, Fs)
% function phi = doa_gcc2(alg,x1,x2,dx,N,Fs)
%
% Estima a direção de chegada (azimuth em graus) para arranjo de mics.
% unidimensional usando a função Correlação Cruzada Generalizada (GCC)
%
% Código adaptado de https://www.nt.tuwien.ac.at/staff-pages/gerhard-doblinger-speaker-tracking-aray/
%
% Referência:
% G. Doblinger:
% "Localization and Tracking of Acoustical Sources";
% in: "Topics in Acoustic Echo and Noise Control", Springer, Berlin - Heidelberg, 2006, ISBN: 3-540-33212-x, 91 - 120.
%
%

% x1,x2     microphone signals
% dx        microphone distance in meters
% N         signal frame length
% Fs        sampling frquency in Hz
% azimuth   azimuth in degrees
% delay     time delay different between mics


%% Reading and defining parameter
x1 = x1(:);
x2 = x2(:);
dx = abs(dx);
Lov = 4;                            % Number of overlapping samples between frames
M = round(N/Lov);                   % Frame hop size
OV = 4;                             % Oversampling factor for interpolation filter
vs = 340;                           % Acoustic waves propagation speed
alpha = 0.7;                        % Forgetting factor of spectral power averaging
doa_threshold = 0.10;               % Minimum correlation value to compute delay

Nx = min(length(x1),length(x2));    % Data vector length
Ndata = 1 + ceil((Nx-N)/M);         % Number of windows processed and thus of calculated delay data points
Nf = N;                             % FFT data points
Nfh = round(Nf/2) + 1;              % FFT data points (half spectrum)
maxDelay = 2 + ceil( dx/vs*Fs );    % Max. delay (delay offset to obtain overall positive
                                    % delays in Cxy)

Nfo = OV * N;                       % Oversampled IFFT length
Ndo = OV * maxDelay;                % Max. delay index for oversampled data
L = 2*maxDelay;                     % Delay range for oversampled data
                                    
                                    
%% Time window function
a = 0.54;
b = -0.46;
w = 2*sqrt(M/N)/(sqrt(4*a^2+2*b^2))*(a + b*cos(pi/N*(2*[0:N-1]'+1))); %Hamming

%% Initializing variables
Cmat = zeros(L*OV, Ndata);
Sxy = zeros(Nfh,1);         % Cross Power Spectrum
delay = zeros(Ndata,1);     % Estimated Time Delay Difference (TDD)
delay_old = Ndo;            % First delay value
m = 0;

%% PHAT-GCC algorithm
for n = 1:M:Nx-N+1  % First window sample index
    m = m+1;        % Window count
    n1 = n:n+N-1;   % Window range
    
    % Compute FFT
    X1 = fft(x1(n1).*w,Nf);
    X2 = fft(x2(n1).*w,Nf);
    X1 = X1(1:Nfh); 
    X2 = X2(1:Nfh);
    
    % Cross-power Spectrum estimation
    Sxy = alpha*Sxy + (1 - alpha)*X1 .* conj(X2);
    
    % IFFT to compute Cross-correlation
    Cxy = OV*real( ifft( Sxy./(abs(Sxy)+1e-4 ), Nfo) ); % ifft zero pads function since Nfo exceeds Sxy length (oversampling)
    Cxy = [Cxy(Nfo-Ndo+1:Nfo) ; Cxy(1:Ndo)];            % Selects correlations within valid delay range
    Cmat(:, m) = Cxy;                                    
    [Cxymax, imax] = max(Cxy); % Finds correlation peak for instant m
    
    % Checks if correlation is valid or noise
    if Cxymax > doa_threshold
       delay(m) = imax-1;
       delay_old = delay(m);
    else
       delay(m) = delay_old;
    end
  
end

t = M/Fs*[0:length(delay)-1];                       % Time vector
tau = linspace(-maxDelay/Fs, maxDelay/Fs, L*OV);    % Delay vector
delay = (delay(1:m+1)-Ndo)/(OV*Fs);                 % Correct delay offset
azimuth = 180/pi*real(acos(vs/dx*delay));           % Azimuth from delay (LINEAR HORIZ. ARRAY)