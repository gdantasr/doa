function y = passby( a, b, speed, x, dist_mic_mic, heigth_mic, sinal)
%   Descri��o matem�tica do �ngulo azimutal ao longo do tempo formado pela
%   posi��o de um ve�culo se deslocando em velocidade constante em uma
%   trajet�ria paralela ao plano do observador (arranjo de microfones).
%   Mais detalhes desse modelo em [1].

%   [1] ROCHA, Gabriela Dantas. Estima��o Das Dire��es De Chegada De Fontes 
%   Sonoras Veiculares Usando Arranjo De Microfones. 2018. TCC (Gradua��o) 
%   - Curso de Engenharia Eletr�nica e de Computa��o, UFRJ, Rio de Janeiro, 2018.

vs = 340;  % Velocidade de propaga��o do som

y = sinal * ( sqrt( ( ( (x - a) .* speed / 3.6) + dist_mic_mic ).^2 + b^2 + heigth_mic^2 ) - sqrt( ( (x - a) .* speed / 3.6 ).^2 + (b)^2 + heigth_mic^2 )) / vs;

end

