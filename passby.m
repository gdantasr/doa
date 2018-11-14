function y = passby( a, b, speed, x, dist_mic_mic, heigth_mic, sinal)
%   Descrição matemática do ângulo azimutal ao longo do tempo formado pela
%   posição de um veículo se deslocando em velocidade constante em uma
%   trajetória paralela ao plano do observador (arranjo de microfones).
%   Mais detalhes desse modelo em [1].

%   [1] ROCHA, Gabriela Dantas. Estimação Das Direções De Chegada De Fontes 
%   Sonoras Veiculares Usando Arranjo De Microfones. 2018. TCC (Graduação) 
%   - Curso de Engenharia Eletrônica e de Computação, UFRJ, Rio de Janeiro, 2018.

vs = 340;  % Velocidade de propagação do som

y = sinal * ( sqrt( ( ( (x - a) .* speed / 3.6) + dist_mic_mic ).^2 + b^2 + heigth_mic^2 ) - sqrt( ( (x - a) .* speed / 3.6 ).^2 + (b)^2 + heigth_mic^2 )) / vs;

end

