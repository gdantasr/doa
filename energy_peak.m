function n90 = energy_peak( signal, fs )
% Retorna n90, o índice correspondente ao instante da passagem do veículo 
% em frente ao arranjo de microfones

timeWindow = 10e-3; % Comprimento da janela da média móvel em seg.
sampWindow = floor( timeWindow * fs ); % Número de amostras da janela
meanSignal = mov_mean(signal, sampWindow); % Média móvel

energy = meanSignal.^2;

[~, n90] = max (energy);

end