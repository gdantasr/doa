function n90 = energy_peak( signal, fs )
% Retorna n90, o �ndice correspondente ao instante da passagem do ve�culo 
% em frente ao arranjo de microfones

timeWindow = 10e-3; % Comprimento da janela da m�dia m�vel em seg.
sampWindow = floor( timeWindow * fs ); % N�mero de amostras da janela
meanSignal = mov_mean(signal, sampWindow); % M�dia m�vel

energy = meanSignal.^2;

[~, n90] = max (energy);

end