function [phi_d, delay, t, tau, Psum] = doa_itd(x1,x2,dx,N,Fs)
%function phi = doa_itd(x1,x2,dx,N,Fs,alpha)
%
% direction estimation (azimuth phi) for 1 dim. microphone array
% using interaural time differences (ITD) map
% resolution set to 2.5 deg. (change parameter dphi and doa_threshold)
%
% (for use in oversampled filterbank)
%
% x1,x2   microphone signals
% dx      microphone distance in meters
% N       signal frame length (512, if omitted)
% Fs      sampling frquency in Hz (16000, if omitted)
% alpha   forgetting factor of ITD maps (0.9 if omitted) 
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
   
% algorithm: C. Liu et al, "Localization of multiple sound sources with two microphones",
%            J.Acoust. Soc. Am. 108 (4), Oct. 2000, pp. 1888-1905.
%
% no filtering of Psum, no threshold applied to ITD map (as suggested by authors)

show_final_map = 0;                          % disable plot of ITD map   

alpha = 0.75;

doa_threshold = 8.5/(1-alpha);               % threshold for speech activity
dphi = 1;%2.5;                                  % azimuth resolution in deg.
Lov = 4;
M = round(N/Lov);                            % frame hop size

% create time window function

a = 0.54;
b = -0.46;
w = 2*sqrt(M/N)/(sqrt(4*a^2+2*b^2))*(a + b*cos(pi/N*(2*[0:N-1]'+1)));

x1 = x1(:);
x2 = x2(:);
Nx = min(length(x1),length(x2));
Ndata = 1+ceil((Nx-N)/M);
Nf = N;
Nfh = Nf/2+1;

phi_d = zeros(Ndata,1);
delay = zeros(Ndata,1);
I = floor(180/dphi);                         % number of azimuths (dphi resolution)
I = I + 1 + mod(I,2);                        % make I odd
imax_old = floor(I/2);

Psum = zeros(Ndata,I-2);
Pmap = zeros(Nfh,I);                         % ITD map
vs = 340;                                    % acoustic waves propagation speed
dx = abs(dx);
Ti = dx/(2*vs)*sin(pi/(I-1)*[0:I-1] - pi/2); % time delay vector
phi = linspace(0,180,I);
f = linspace(0,Fs/2,Nfh)';
e1 = exp(-j*2*pi*f*Ti);                      % matrix of phase factors (first signal)
e2 = exp(-j*2*pi*f*Ti(end:-1:1));            % vector of phase factors (second signal)
t = M/Fs*[0:Ndata-1];

m = 0;
for n = 1:M:Nx-N+1
   m = m+1;                                  % frame counter
   n1 = n:n+N-1; 
   X1 = fft(x1(n1).*w,Nf);
   X2 = fft(x2(n1).*w,Nf);
   X1 = X1(1:Nfh)*ones(1,I);                 % use half of the spectra (real-valued signals)
   X2 = X2(1:Nfh)*ones(1,I);
   
   X1 = X1 .* e1;                            % matrix of "delayed" spectrum
   X2 = X2 .* e2;
   dX = abs(X1-X2);                          % coincidence detection cost function eq. (6) 
   [dum,imin] = min(dX.');                   % eq. (5) 
   Pmap = alpha*Pmap;                        % exponentially weighted (over time) histogram
   for k = 1:Nfh
      it = imin(k);
      Pmap(k,it) = Pmap(k,it) + 1;           % update map of histogram
   end
   Pmap(find(Pmap<1)) = 0;                   % this may eliminate some amount of spurious data
   Psum(m,:) = sum(Pmap(:,2:I-1));           % sum-up (average over frequency)   
                                             % eliminate minima detection erros at phi = 0, 180°
                                             % (at frequencies with no spectral components)
   [Pmax,imax] = max(Psum(m,:));             % locate maximum
   if Pmax > doa_threshold
      phi_d(m) = phi(imax);
      delay(m) = 2*Ti(imax);
      imax_old = imax;
   else
      phi_d(m) = phi(imax_old);
      delay(m) = 2*Ti(imax_old);
   end
end

t = M/Fs*[0:m-1];
delay = 2*delay;
maxDelay = 2 + ceil( dx/vs*Fs );
tau = linspace(-maxDelay/Fs, maxDelay/Fs, length(Ti));    % Delay vector
tau = tau(2:end-1);
Psum = flipud(Psum.');

% % plot phi-trace
% 
% close all
% pos = [0.01 0.5 0.49 0.42];
% figurebackcolor = 'black';
% fp1 = figure('numbertitle','off','name','DOA estimation',...
%              'Units','normal','Position',pos);
% colordef(fp1,figurebackcolor);
% t = M/Fs*[0:m-1];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ind2,phi2] = TheoPhi_and_Espectrogram(v,t90,x1,Fs);
% plot(t,phi_d(1:m),ind2,phi2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% xlabel('Time t in sec.'), ylabel('Azimuth \phi in deg.');
% title(['Method: ITD, ','d = ', num2str(100*dx),' cm,  ', 'F_s = ', num2str(Fs), ' Hz', ', V = ', num2str(v), ' Km/h']);
% grid on
% axis tight
% 
% % plot final Pmap
% show_final_map=1;
% if show_final_map
%    pos = [0.01 0.05 0.49 0.38];
%    fp1 = figure('numbertitle','off','name','ITD map',...
%                 'Units','normal','Position',pos);
%    colordef(fp1,figurebackcolor);
% 
%    imagesc(phi,f,Pmap);
%    set(gca,'YDir','normal');
%    ylabel('Freq. in Hz');
%    xlabel('Azimuth \phi in deg.');
%    title(sprintf('ITD histogram map at time = %3.2f sec', t(end)));
% end
% 
% plot Psum over time

% figurebackcolor = 'black';
% pos = [0.505 0.5 0.49 0.42];
% fp1 = figure('numbertitle','off','name','ITD map',...
%              'Units','normal','Position',pos);
% colordef(fp1,figurebackcolor);
% 
% Psum = Psum(1:m,:)';
% imagesc(t,phi,Psum), colormap(cool), colorbar
% set(gca,'YDir','normal');
% xlabel('Time t in sec.');
% ylabel('Azimuth \phi in deg.');
% title('Frequency averaged ITD histograms vs. time');
% pause


