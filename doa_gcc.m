function [azimuth, delay, t, tau, Cmat] = doa_gcc(x1, x2, dx, N, Fs)
%function phi = doa_gcc2(alg,x1,x2,dx,N,Fs)
%
% direction estimation (azimuth phi) for 1 dim. microphone arrays
% using generalized cross correlation and speech activity detection
%

% x1,x2     microphone signals
% dx        microphone distance in meters
% N         signal frame length
% Fs        sampling frquency in Hz
% azimuth   azimuth in degrees

Lov = 4;    % Number of overlapping samples between frames
M = round(N/Lov);   % frame hop size

% Time window function
a = 0.54;
b = -0.46;
w = 2*sqrt(M/N)/(sqrt(4*a^2+2*b^2))*(a + b*cos(pi/N*(2*[0:N-1]'+1))); %Hamming

x1 = x1(:);
x2 = x2(:);
Nx = min(length(x1),length(x2)); % Data vector length
Ndata = 1+ceil((Nx-N)/M); % Number of windows processed and thus of calculated delay data points
Nf = N;

% Initializing vectors
Nfh = round(Nf/2) + 1;  % FFT data points
Sxy = zeros(Nfh,1);     % Cross Power Spectrum
delay = zeros(Ndata,1); % Estimated Time Delay Difference (TDD)

dx = abs(dx);
OV = 4;                 % Oversampling factor for interpolation filter
vs = 340;               % Acoustic waves propagation speed
Nd = 2+ceil(dx/vs*Fs);  % Max. delay (delay offset to obtain overall positive
                        % delays in Cxy)
Nfo = OV*Nf;
Ndo = OV*Nd;

L = 2*Nd;

alpha = 0.7;    % forgetting factor of spectral power averaging
delay_old = Ndo;
doa_threshold = 0.10;

Cmat = zeros(L*OV,Ndata);

m = 0;
% PHAT-GCC algorithm
for n = 1:M:Nx-N+1
    m = m+1;
    n1 = n:n+N-1; 
    X1 = fft(x1(n1).*w,Nf);
    X2 = fft(x2(n1).*w,Nf);
    
    X1 = X1(1:Nfh); 
    X2 = X2(1:Nfh);
    Sxy = alpha*Sxy + (1 - alpha)*X1 .* conj(X2);
    Cxy = OV*real( ifft( Sxy./(abs(Sxy)+1e-4 ), Nfo) ); % ifft zero pads function since Nfo exceeds Sxy length (oversampling)
    Cxy = [Cxy(Nfo-Ndo+1:Nfo) ; Cxy(1:Ndo)];
    Cmat(:,m) = Cxy; 
    [Cxymax,imax] = max(Cxy);
    
    if Cxymax > doa_threshold
       delay(m) = imax-1;
       delay_old = delay(m);
    else
       delay(m) = delay_old;
    end
  
end

t = M/Fs*[0:length(delay)-1];
tau = linspace(-Nd/Fs,Nd/Fs,L*OV);
delay = (delay(1:m+1)-Ndo)/(OV*Fs);            % correct delay offset
azimuth = 180/pi*real(acos(vs/dx*delay));