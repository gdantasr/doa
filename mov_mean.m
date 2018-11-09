function meanSignal = mov_mean( signal, wLen )
%   Retorna m�dia do sinal calculada em uma janela deslizante w de wLength
%   colunas, centrada na posi��o atual se o comprimento
%   for �mpar. Se wLen for par, aumenta em 1 o comprimento
%   A m�dia � calculada


% Resto da divis�o por 2 para checar se � �mpar
if mod (wLen, 2) == 0
    wLen = wLen + 1; % Se for par soma 1 ao comprimento
end

signalLen = length(signal);
meanSignal = zeros(size(signal)); % Inicializa vetor

for k = 1:length(signal)
    
    % �ndices dos limites da janela
    kDown = max ( [1, k - floor(wLen)] ); % �ndice inferior � no m�nimo 1
    kUp =  min ( [signalLen, k + floor(wLen)] ); % �ndice superior � no m�ximo o comprimento do sinal
    
    % Calcula m�dia
    meanSignal(k, :) = mean(signal(kDown : kUp, :));
   
end


end