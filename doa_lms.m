function [phi, delay, t, tau, Wmat] = doa_lms(x1,x2,dx,N,M,mu,Fs)

%[phi, tau, Wmat] = doa_fastlmsfullBand(x1,x2,dx,N,M,mu,Fs)
%function phi = doa_fastlms(x1,x2,dx,N,M,mu,Fs)
%
% direction estimation (azimuth phi) for 1 dim. microphone arrays
% using adaptive filter (fast block LMS algorithm with frequency dep. mu)
%
% x1,x2   microphone signals
% dx      microphone distance in meters
% N       adaptive filter length (FIR filter)
% M       delay is computed every M samples (M must be a multiple of N)  
% mu      step size of adaptive algorithm (typ < 0.5)
% Fs      sampling frquency in Hz (16000, if omitted)
% phi     azimuth in degrees
%
%   Copyright 2006 Gerhard Doblinger, Vienna University of Technology
%   g.doblinger@tuwien.ac.at
%   http://www.nt.tuwien.ac.at/about-us/staff/gerhard-doblinger/
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

% ref.: F.A.Reed, P.L.Feintuch, N.J.Bershad, "Time delay estimation using the LMS adaptive
%       filter - static behavior", IEEE Trans. Acoust., Speech, Signal Proc., vol. ASSP-29,
%       no. 3, June 1981, pp 561-571.

vs = 340;                 % acoustic waves propagation speed
dx = abs(dx);
Nd = 2+ceil(dx/vs*Fs);    % max. delay between mics in samples
Lov = 4;                  % oversampling factor (for u2 vector)
Fs1 = Lov*Fs;
Ndo = Nd*Lov;

doa_threshold = 0.1;      % speech activity threshold
                          % CHANGE, if necessary

x1 = x1(:);
x2 = x2(:);
Nx = min(length(x1),length(x2));

% set M to a multiple of filter lenght N (required by block processing)

Mb = max(1,round(M/N));
M = Mb*N;
Nf = 2*N;                 % FFT length

if Nx <= Nf
   disp('INFO: FFT length must be less than signal length');
   return
end

% init. variables

Nh = floor(N/2);          % delay of first microphone signal
W = zeros(Nf,1);          % weight vector in frequency domain
Nb = ceil((Nx-Nf)/N);
delay = zeros(Nb,1);
delay_old = 0;
Wmat = zeros(2*Ndo,Nb);

Sx = zeros(Nf,1);         % spectral power used with step size mu
alpha = 0.2;              % forgetting factor of spectral power averaging
alpha1 = 1-alpha;

% adaptive filter loop

k = 0;
mb = 0;
for n = 1:N:Nx-Nf
   mb = mb+1;
   X2 = fft(x2(n:n+Nf-1));
   y = real(ifft(W.*X2));
   e = x1(n+Nh+1:n+Nh+N) - y(N+1:Nf);
   E = fft([zeros(N,1) ; e]);
   Sx = alpha*Sx + alpha1*abs(X2).^2;
   W = W + (mu ./ (Sx + eps)) .* conj(X2) .* E; 
   if mod(mb,Mb) == 0          % determine delay every Mb blocks only
      w = real(ifft(W));
      w1 = w(Nh-Nd:Nh+Nd-1);   % weight vector to be used to find delay
      w1 = resample(w1,Lov,1); % interpolate
      k = k+1;
      Wmat(:,k) = w1;          % filter coefficient map
      [wmax,dmax] = max(w1);
      del = dmax-Ndo;
      if wmax < doa_threshold  
         delay(k) = delay_old; % use old delay during speech pauses 
      else
         delay(k) = del;
         delay_old = del;
      end
   end
end

delay = delay(1:k);
Wmat = Wmat(:,1:k);
tau = linspace(-Nd/Fs,Nd/Fs,2*Ndo);
phi = 180/pi*real(acos(vs/dx*delay/Fs1));
t = linspace(0,(Nx-Nf-1)/Fs,k);

%disp(sprintf('final delay = %2.4e msec (%d samples)', ...
%             1000*delay(end)/Fs1, delay(end)));
%phi = 180/pi*real(acos(vs/dx*delay/Fs1));
%disp(sprintf('final azimuth = %2.4e deg', phi(end)));


% % plot results (phi, final weight vector)
% 
% close all
% figurebackcolor = 'black';
% pos = [0.01 0.5 0.49 0.42];
% fp1 = figure('numbertitle','off','name','DOA estimation',...
%              'Units','normal','Position',pos);
% colordef(fp1,figurebackcolor);
% t = linspace(0,(Nx-Nf-1)/Fs,k);
% 
% grid on, xlabel('Time t in sec.'), ylabel('Azimuth \Phi in deg.');
% % title(['Method: FLMS, ','d = ', num2str(100*dx),' cm,  ', 'F_s = ', num2str(Fs), ...
% %        ' Hz, N = ',int2str(N), ', M = ', int2str(M), ', V = ', num2str(v), ' Km/h']);
% grid on;
% axis tight;
% 
% pos = [0.5055 0.5 0.49 0.42];
% fp1 = figure('numbertitle','off','name','DOA estimation',...
%              'Units','normal','Position',pos);
% colordef(fp1,figurebackcolor);
% tau = 1000*linspace(-Nd/Fs,Nd/Fs,2*Ndo);
% imagesc(t,tau,Wmat), colormap cool, colorbar;
% set(gca,'Ydir','normal');
% xlabel('Time t in sec.'), ylabel('Delay \tau in msec.');
% title('Adaptive filter coefficients w[t,\tau]');
% 
% pos = [0.01 0.06 0.49 0.36];
% fp1 = figure('numbertitle','off','name','DOA estimation',...
% 	     'Units','normal','Position',pos);
% colordef(fp1,figurebackcolor);
% 
% plot(tau,w1,tau(dmax),wmax,'or'), grid on;
% axis tight;
% xlabel('\tau in msec.'), ylabel('w[\tau]');
% title(['final weight vector, oversampling factor = ', ...
%        int2str(Lov)]);


