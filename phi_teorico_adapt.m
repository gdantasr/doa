function [phi, phi2, dt, dt2] = phi_theorico (v, t90_inf, t90_sup, sentido, dist_source_mic, dist_source_source, heigth_mic, T)

% phi e phi2 correspondem as curvas teoricas da DOA para as duas fontes
% phi3 desconsidera a presenca das multiplas fontes, usando o ponto central
% entre elas

sound_speed = 340;

t=0;
ind2=0:0.01:T;
for t2=0:0.01:T
    t=t+1;
    l_inf(t) = ((t2 - t90_inf)*v/3.6); % Distancia em METROS.
    l_sup(t) = ((t2 - t90_sup)*v/3.6);
    d = 0.2;
    d1 = sqrt( (l_sup(t))^2 + (heigth_mic)^2 + (dist_source_mic)^2 );
    d2 = sqrt( (l_sup(t) + d)^2 + (heigth_mic)^2 + (dist_source_mic)^2 );
    delta(t) = d2 - d1;
    d3 = sqrt( (l_inf(t))^2 + (heigth_mic)^2 + (dist_source_mic)^2 );
    d4 = sqrt( (l_inf(t) + d)^2 + (heigth_mic)^2+(dist_source_mic)^2 );
    delta2(t) = d4 - d3;
    %d5 = sqrt( (l(t))^2 + (heigth_mic)^2 + (dist_source_mic)^2 );
    %d6 = sqrt( (l(t) + d)^2 + (heigth_mic)^2+(dist_source_mic)^2 );
    %delta3(t) = d6 - d5;
    
    if sentido == '0° -> 180°'
        delta(t) = -delta(t);
        delta2(t) = -delta2(t);
        %delta3(t) = -delta3(t);
    end
     
    dt(t) = delta(t)/sound_speed;
    dt2(t) = delta2(t)/sound_speed;
    %dt3(t) = delta3(t)/sound_speed;
    phi(t) = acos( delta(t)/d ) * 180/pi;
    phi2(t) = acos( delta2(t)/d ) * 180/pi;
    %phi3(t) = acos( delta3(t)/d ) * 180/pi;
end

end

