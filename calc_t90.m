function n90 = calc_t90(x1,fs)
%function [t90, Nspec90, N_spec_win, f_spec, media_s_spec] = calc_t90(x1,fs)

N_janela = 1024;
[S, F, t, P] = spectrogram(x1,hamming(N_janela),N_janela/2,N_janela,fs,'yaxis');

totalBins = N_janela/2;
f_min = 500 * totalBins / fs * 2 + 1;
f_max = 1000 * totalBins / fs * 2 + 1;

% Media da energia para as frequencias que sofrem grande influencia do
% ruido veicular. Restricao das freqs para evitar outros ruidos.
P_avg = mean(P(f_min : f_max, :));

[~, maxEnergyIndex] = max(P_avg);   % Indice de maior energia da FFT

t90 = (maxEnergyIndex * t(end)) / length(t);
n90 = round(t90*fs);
