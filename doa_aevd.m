function [phi, delay, t, tau, U] = doa_aevd(x1, x2, dx, N, M, mu, Fs)
%function phi = doa_aevd2(x1,x2,dx,N,M,mu,Fs)
%
% direction estimation (azimuth phi) for 1 dim. microphone arrays
% using adaptive eigenvalue decomposition
% version with re-initialization every M samples (to ensure tracking)
% and fast block LMS algorithm
%
% x1,x2   microphone signals
% dx      microphone distance in meters
% N       eigenfilter length (FIR filter)
% M       re-initialization period in samples (must be a multiple of N)  
% mu      step size of adaptive algorithm
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

% ref.: J.Benesty, "Adaptive eigenvalue decomposition algorithm for passive
%       acoustic source localization", J.Acoust.Soc.Am. 107, Jan. 2000,
%       pp. 384-391.

vs = 340;                 % acoustic waves propagation speed
dx = abs(dx);
Nd = 2+ceil(dx/vs*Fs);    % max. delay between mics in samples
Lov = 4;                  % oversampling factor (for u2 vector)
Fs1 = Lov*Fs;
Ndo = Nd*Lov;

% doa_threshold = -0.09;    % speech activity threshold
%                           % CHANGE, if necessary
doa_threshold = -0.009;    % speech activity threshold


x1 = x1(:);
x2 = x2(:);
Nx = min(length(x1),length(x2));

% set M to a multiple of filter lenght N (required by block processing)

Mb = max(1,round(M/N));
M = Mb*N;
Nf = 2*N;                 % FFT length   

% init. eigenvector u

Nh = floor(N/2)+1;
u0_1 = zeros(N,1);
u0_1(Nh) = 1;
U0_1 = fft(u0_1,Nf);
U0_2 = zeros(Nf,1);
U1 = U0_1;
U2 = U0_2;
Nb = ceil((Nx-N)/M);
delay = zeros(Nb,1);
delay_old = 0;
U = zeros(2*Ndo,Nb);

Sx1 = zeros(Nf,1);        % spectral power used with step size mu
Sx2 = zeros(Nf,1);
alpha = 0.2;              % forgetting factor of spectral power averaging  
alpha1 = 1-alpha;

% loop to compute eigenvector u corresponging to zero eigenvalue

k = 0;
mb = 0;
for n = 1:N:Nx-Nf
   mb = mb+1;
   m = n:n+Nf-1;
   X1 = fft(x1(m));
   X2 = fft(x2(m));
   e = real(ifft(U1.*X1+U2.*X2));
   E = fft([zeros(N,1) ; e(N+1:Nf)]);
   Sx1 = alpha*Sx1 + alpha1*abs(X1).^2;
   Sx2 = alpha*Sx2 + alpha1*abs(X2).^2;
   U1 = U1 - (mu ./ (Sx1+eps)) .* conj(X1) .* E;
   U2 = U2 - (mu ./ (Sx2+eps)) .* conj(X2) .* E;
   if mod(mb,Mb) == 0           % find delay, and restart adaptive filter
      u2 = real(ifft(U2));      % eigenvector to be used to find delay
      u2 = u2(Nh-Nd:Nh+Nd-1);
      u2 = resample(u2,Lov,1);  % interpolate
      k = k+1;
      U(:,k) = u2;
      [umin,dmin] = min(u2);
      del = dmin-Ndo;
      if umin > doa_threshold   % signal to weak (e.g. speech pauses)
         delay(k) = delay_old;
      else
         delay(k) = del;
         delay_old = del;
      end
      U1 = U0_1;
      U2 = U0_2;
   end
end

tau = linspace(-Nd/Fs,Nd/Fs,2*Ndo);
delay = delay(1:k);
phi = 180/pi*real(acos(vs/dx*delay/Fs1));
t = linspace(0,(Nx-N-1)/Fs,length(phi));

disp(sprintf('final delay = %2.4e msec (%d samples)', ...
             1000*delay(end)/Fs1, delay(end)));
disp(sprintf('final azimuth = %2.4e deg', phi(end)));

% plot results (phi, final eigenvector)

% close all
% figurebackcolor = 'black';
% pos = [0.01 0.5 0.49 0.42];
% fp1 = figure('numbertitle','off','name','DOA estimation',...
%              'Units','normal','Position',pos);
% colordef(fp1,figurebackcolor);
% t = linspace(0,(Nx-N-1)/Fs,length(phi));
% grid on, xlabel('Time t in sec.'), ylabel('Azimuth \Phi in deg.');
% title(['Method: AEVD, ','d = ', num2str(100*dx),' cm,  ', 'F_s = ', num2str(Fs), ...
%        ' Hz, N = ',int2str(N), ', M = ', int2str(M), ', V = ', num2str(v), ' Km/h']);
% grid on;
% axis tight;
% 
% pos = [0.5055 0.5 0.49 0.42];
% fp1 = figure('numbertitle','off','name','DOA estimation',...
%              'Units','normal','Position',pos);
% colordef(fp1,figurebackcolor);
% tau = 1000*linspace(-Nd/Fs,Nd/Fs,2*Ndo);
% imagesc(t,tau,-U), colormap cool, colorbar;
% set(gca,'Ydir','normal');
% axis tight;
% xlabel('Time t in sec.'), ylabel('Delay \tau in msec.');
% title('Impulse response h_1[t,\tau]');
% 
% pos = [0.01 0.06 0.49 0.36];
% fp1 = figure('numbertitle','off','name','DOA estimation',...
%              'Units','normal','Position',pos);
% colordef(fp1,figurebackcolor);
% 
% plot(tau,u2,tau(dmin),umin,'or'), grid on;
% axis tight;
% xlabel('\tau in msec.'), ylabel('u_2[\tau]');
% title(['final eigenvector u_2, oversampling factor = ', ...
%        int2str(Lov)]);


