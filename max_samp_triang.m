function [ num_amostras ] = max_samp_triang( angulo, velocidade,fs )
% Função que dado um angulo correspondente a um deslocamento do veículo, calcula o número
% máximo de amostras compreendidas nesse deslocamento.
% 
% angulo : Medida em graus do intervalo onde a direção de chegada será
% considerada constante
% velocidade : Velocidade do veículo medida em km/h

distancia = 3.6; % Distancia medida em metros entre o veículo e o array de microfones 
%fs = 44100; % Frequência de amostragem em Hz

velocidade = velocidade/3.6; % Velocidade em m/s
deslocamento = 2*distancia*sind(angulo/2)/cosd(angulo/2); 
t_desloc = (deslocamento/velocidade); % Tempo em segundos decorrido durante o deslocamento
num_amostras = floor(t_desloc*fs);

end