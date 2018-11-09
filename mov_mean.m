function meanSignal = mov_mean( signal, wLen )
%   Retorna média do sinal calculada em uma janela deslizante w de wLength
%   colunas, centrada na posição atual se o comprimento
%   for ímpar. Se wLen for par, aumenta em 1 o comprimento
%   A média é calculada


% Resto da divisão por 2 para checar se é ímpar
if mod (wLen, 2) == 0
    wLen = wLen + 1; % Se for par soma 1 ao comprimento
end

signalLen = length(signal);
meanSignal = zeros(size(signal)); % Inicializa vetor

for k = 1:length(signal)
    
    % Índices dos limites da janela
    kDown = max ( [1, k - floor(wLen)] ); % Índice inferior é no mínimo 1
    kUp =  min ( [signalLen, k + floor(wLen)] ); % Índice superior é no máximo o comprimento do sinal
    
    % Calcula média
    meanSignal(k, :) = mean(signal(kDown : kUp, :));
   
end


end