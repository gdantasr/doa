function [y, RdB] = adapt_beam(Xin,mics,phi_doa,mu,Gw,f,Fs)
%function y = adaptive_beam(X,mics,phi_look,mu,Gw,f,Fs)
%
% adaptive beamformer using an FFT filter bank with/without DOA estimation and tracking
% (frequency domain version of Frost beamformer)
% version to plot array pattern (image and 3d plot)
%
% X          matrix of microphone signals (one column for each channel)
% mics       microphone placement
%            (first col = x, second col = y coord.)
% phi_look   azimuth in deg. of array look direction (0 ... 180°)
%            if phi_look=[] (or omitted), then automatic speaker tracking is used 
% mu         step size of adaptive algorithm (default 0.01)
% Gw         Gw parameter in dB limit at low frequencies (default 4)
% f          frequency of array pattern in Hz (default 1000 Hz) 
% Fs         sampling frequency in Hz (default 16000 Hz)
% y          beamformer output signal
%
% example:   y = adaptive_beam(X,mics);    % automatic speaker tracking
%            y = adaptive_beam(X,mics,90); % look direction = 90°    

if nargin == 0
   help adaptive_beam
   return
elseif nargin == 2
   phi_look = [];
   mu = 0.01;
   Gw = 4;
   f = 500;
   Fs = 8000;
elseif nargin == 3
   mu = 0.01;
   Gw = 4;
   f = 500;
   Fs = 8000;
elseif nargin == 4
   Gw = 4;
   f = 500;
   Fs = 8000;
elseif nargin == 5
   f = 500;
   Fs = 8000;
elseif nargin == 6
   Fs = 8000;
end

if f >= Fs/2
   error('f >= Fs/2');
end

[Nx,N] = size(Xin);
[Nc,dum] = size(mics);
if Nc ~= N
   error('number of mics does not match number of input signals');
end

Nfft = 64;                          % FFT length = time window length
% Nfh = Nfft/2+1;
M = Nfft/4;                          % frame hop size 
disp(sprintf('Nmics = %d, Nfft = %d, Fs = %5.0f Hz',N, Nfft, Fs));

% set max. input signal frequency according to min. array spacing

vs = 340;                            % acoustic waves propagation speed
dmin = min(norm(diff(mics(:,1))));
fmax = Fs/2;
Nfmax = round(Nfft*fmax/Fs);
disp(sprintf('Nfmax = %d, fmax = %5.0f Hz', Nfmax, fmax));

% compute parameters of starting solution (at azimuth phi_d)  

phi_look = abs(phi_doa(1));
theta_d = 90*ones(1, length(phi_look))
phi_d = phi_look;
[b,Wc,P] = compute_params(mics,vs,Nfft,Nfmax,theta_d,phi_d,Gw,Fs);

Ndata = 1+ceil((Nx-Nfft)/M);
%phi_doa = phi_look(1)*ones(Ndata,1);
disp(sprintf('no speaker tracking, array look direction = %d deg.', phi_look));

% variabels of array pattern

phi = pi/180*[0:180];
Nphi = length(phi);
er = mics*sin(pi/180*theta_d(1))*[cos(phi) ; sin(phi)];
mf1 = round(Nfft*f/Fs);              % nearest frequency index
beta = 2*pi*mf1*Fs/(vs*Nfft);
D1 = exp(j*beta*er);                 % matrix of steering vectors
Ndata = 1+ceil((Nx-Nfft)/M);
RdB = zeros(Ndata,Nphi);

% beam look direction (from DOA estimation)
phi_doa = resample(phi_doa, Ndata-1, length(phi_doa));

% filter bank 

W = Wc;                              % initial adaptive filter coefficients
V = zeros(size(Wc));
y = zeros(Nx,1);
h = hanning(Nfft);
H = h(:) * ones(1,N);                % matrix of time windows

disp('computing beam pattern, and adaptive array output...');
drawnow

phi_old = phi_doa(1); % Angulo inicial
n = 0;
for  m = 1:M:Nx-Nfft+1
    
    % compute filter bank output signal
    n = n+1;
    m1 = m:m+Nfft-1; 
    X = fft(Xin(m1,:) .* H).';        % apply FFT to weighted frames of all mics
    Y = sum(conj(W) .* X(:,1:Nfmax)); % apply adaptive filter coefficients, sum up FFT bins
    yb = real(ifft(Y,Nfft));
    y(m1) = y(m1) + yb(:);            % overlap-add
   
    % compute new starting solution for new estimated direction
    if abs(phi_doa(n)-phi_old) >= 2           % new start with new parameters  
        [b,Wc,P] = compute_params(mics,vs,Nfft,Nfmax,theta_d,phi_doa(n),Gw,Fs);
        phi_old = phi_doa(n);    % save DOA estimate
    end
 
    % update coefficient matrix W (constrained LMS algorithm)
   
    for k = 1:Nfmax
       V(:,k) = P(:,:,k) * (V(:,k) - mu/(norm(X(:,k))^2+eps)*X(:,k)*conj(Y(k)));
       no = norm(V(:,k));
       if no > b(k)                    % norm constraint
           V(:,k) = V(:,k)*b(k)/no;
       end
    end
    W = Wc + V;
   
    % compute array pattern
    R = abs(W(:,mf1)'*D1).^2;
    RdB(n,:) = R/max(R);

end

y = y/Nfft*M*4;
% audio = audioplayer(y,Fs);
% playblocking(audio)
audioname = strcat('beam_fs',num2str(Fs), '_f', num2str(f),'.wav');
audiowrite(audioname,y,Fs);



% plot 2-dim. array pattern

close all
pos = [0.01 0.5 0.49 0.4];
figurebackcolor = 'black';
fabf = figure('numbertitle','off','name','Adaptive beamformer I',...
	     'Units','normal','Position',pos, 'DoubleBuffer', 'on');
colordef(fabf,figurebackcolor);
RdB = RdB(1:n,:);
RdB = RdB';
RdB = max(-40,10*log10(RdB+eps));
t = M/Fs*[0:n-1];
map = colormap('jet');
%colormap(flipud(map));

hold on
grid on
imagesc(t, 180/pi*phi, RdB);
colorbar
plot(t, phi_doa, 'y');
hold off
set(gca,'YDir','normal');
ylabel('Azimuth \Phi in deg.');
xlabel('Time t in sec.');
title(['Array pattern in dB at f = ', num2str(f),' Hz, Fs =', num2str(Fs),' Hz, \mu = ', num2str(mu)]);
text(1.04*t(end), 185,'  dB');
if isempty(phi_look)
   text(3, 96, 'Estimated Azimuth', 'Color', 'yellow');
else
   text(3, phi_look(1)+6, 'Look Direction', 'Color', 'yellow');
end
axis([0, t(end),0,180]);

fig = gcf;
fig.InvertHardcopy = 'off';        
filename = strcat('beam_pattern_fs',num2str(Fs), '_f', num2str(f),'.png');
%print(gcf, filename,'-dpng','-r300'); % Salvando imagem em arquivo png
% pause

% 
% figure
% mesh(RdB); hold on
% plot(phi_doa, 'b');


% plot 3-dim. array pattern

pos(1) = pos(1)+0.5;
   fp = figure('numbertitle','off','name','Adaptive beamformer II',...
	       'Units','normal','Position',pos);
colordef(fp,figurebackcolor);

id = 1:50:length(t);
colormap('gray');
surfc(t,180/pi*phi,RdB,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
surfl(t(id),180/pi*phi,RdB(:,id));
shading interp
ylabel('Azimuth \Phi in deg.');
xlabel('Time t in sec.');
zlabel('dB');
title(['Array pattern in dB at f = ', num2str(f),' Hz, \mu = ', num2str(mu)]);
axis([0,t(end),0,180,-40,0]);
view([30,45]);
pause

%---------------------------------------------------------------------------

function w = delaysum_beam(mics,vs,theta_d,phi_d,f);
% function w = delaysum_beam(mics,theta_d,phi_d,f);
%
% compute weights of delay&sum beamformer
%
% w              weight vector of length N (number of mics)
% mics           [xi,yi] microphone cooerdinates (N x 2 matrix)
% vs             speed of spound im m/sec  
% theta_d        desired elevation in deg.
% phi_d          desired azimuth in deg.
% f              frequency in Hz
%
% G.Doblinger, TU-Wien 10-04

if nargin ~= 5
   help delaysum_beam;
   return;
end

theta_d = theta_d*pi/180;
phi_d = phi_d*pi/180;

beta = 2*pi*f/vs;         % wave number
[N,K] = size(mics);
if (K < 2) | (N < 1)
   error('bad matrix of microphine coordinates');
end

if K == 2
   rn = [mics zeros(N,1)];
else
   rn = mics;
end

% steering vector of desired direction

ed = [sin(theta_d)*cos(phi_d); sin(theta_d)*sin(phi_d); cos(theta_d)];
d0 = exp(j * beta * rn * ed);
w = d0 / N;

%---------------------------------------------------------------------------

function [b,Wc,P] = compute_params(mics,vs,Nfft,Nfmax,theta_d,phi_d,Gw,Fs)

% compute parameters of starting solution

[N,dum] = size(mics);                % N = number of microphones
if Nfmax >  Nfft/2+1
   error('Nfmax  > Nfft/2+1');
end
Nco = length(theta_d);               % number of constraints
P = zeros(N,N,Nfmax);
Wc = zeros(N,Nfmax);
b = zeros(Nfmax,1);
if Nco == 1
   for l = 1:Nfmax
      if ((l-1)/Nfft*Fs > 250) & (Gw < -8)  Gw = -8; end
      if ((l-1)/Nfft*Fs > 450) & (Gw < -2)  Gw =  -2; end
      if ((l-1)/Nfft*Fs > 700) & (Gw < 2)  Gw =  2; end
      if ((l-1)/Nfft*Fs > 1000) & (Gw < 4) Gw =  4; end
      if ((l-1)/Nfft*Fs > 2000) & (Gw < 6) Gw =  6; end
      c = delaysum_beam(mics,vs,theta_d,phi_d,(l-1)/Nfft*Fs);
      P(:,:,l) = eye(N) - N*c*c';
      Wc(:,l) = c;
      b(l) = sqrt(10^(-Gw/10) - Wc(:,l)' * Wc(:,l));     
   end
else
   g = [1;zeros(Nco-1,1)];           % vector of constraints
   C = zeros(N,Nco,Nfmax);
   Creg = 1e-10*eye(Nco);            % add regularization term to avoid singular matrix
   for l = 1:Nfmax
      if ((l-1)/Nfft*Fs > 250) & (Gw < -8)  Gw = -8; end
      if ((l-1)/Nfft*Fs > 450) & (Gw < -2)  Gw =  -2; end
      if ((l-1)/Nfft*Fs > 700) & (Gw < 2)  Gw =  2; end
      if ((l-1)/Nfft*Fs > 1000) & (Gw < 4) Gw =  4; end
      if ((l-1)/Nfft*Fs > 2000) & (Gw < 6) Gw =  6; end
      for co = 1:Nco
          C(:,co,l) = N*delaysum_beam(mics,vs,theta_d(co),phi_d(co),(l-1)/Nfft*Fs);
      end
      C1 = C(:,:,l) * (Creg + C(:,:,l)' * C(:,:,l))^-1; 
      P(:,:,l) = eye(N) - C1*C(:,:,l)';
      Wc(:,l) = C1 * g;
      b(l) = sqrt(10^(-Gw/10) - Wc(:,l)' * Wc(:,l));     
   end
end
b = real(b);                         % remove small imaginary part due to rounding

%---------------------------------------------------------------------------




