function [t90, Nspec90, f_spec] = calc_t90(x1,fs)
%function [t90, Nspec90, N_spec_win, f_spec, media_s_spec] = calc_t90(x1,fs)

N_janela = 1024;

[s_spec,f_spec,t_spec] = spectrogram(x1,hamming(N_janela),N_janela/2,N_janela,fs,'yaxis');

for i = 1:length(t_spec)
    max_s_spec(i) = max(s_spec(:,i));
end


N_ov = 4;                       % overlap
N_win = fs;                     % tamanho da janela usada para o intervalo de maxima energia (passagem do carro)
N_hop = round(N_win/N_ov);      % salto da janela
N_samp = length(x1);            % numero total de amostras do sinal reamostrado
tempoTotal = N_samp/fs;

%calcula energia em cada janela e salva num vetor
m=1;
for n = 1:N_hop:N_samp-N_win+1                          
    win_p(m) = sqrt(sumsqr(x1(n:n+N_win-1)));  
    samp_lim(1,m) = n;                                     %amostras que delimitam os intervalos
    samp_lim(2,m) = n+N_win-1;                             %de cada janela no sinal reamostrado   
    m=m+1;
end


[pks,locs] = findpeaks(win_p);  %obtem os picos do vetor de energia
N_picos = length(locs);
limiar = 10;                    %limiar para eliminar ou repor picos 
win_p_aux = win_p;

%ajusta um limiar para deixar passar apenas os picos importantes
qtdeCarros = 1; %%%%%%%%%%%%
while N_picos ~= qtdeCarros 
    win_p_aux = win_p;
    if N_picos < qtdeCarros
        limiar = limiar - 1;
        win_p_aux(win_p_aux<limiar)=0;
    elseif N_picos > qtdeCarros 
        limiar = limiar + 1;
        win_p_aux(win_p_aux<limiar)=0;
    end
    [pks,locs] = findpeaks(win_p_aux);
    N_picos = length(locs);
end


%calcula as amostras do espectrograma que correspondem ao do sinal reamostrado para cada carro a cada iteracao
% Tinha esse for aqui quando eram mais carros
%for i=1:qtdeCarros
    
    Nispec = round(samp_lim(1,locs)*length(max_s_spec)/N_samp);
    Nfspec = round(samp_lim(2,locs)*length(max_s_spec)/N_samp);
    Nspec90_ref = (Nispec+Nfspec)/2;    %referencial para amostra em 90 graus
    Nspec90 = Nfspec;                %condicao inicial para amostra em 90 graus

    %tentativa de aproximar o pico do valor esperado
    erro_spec = 4;
    while abs(Nspec90-Nspec90_ref) > erro_spec
        [max_s_specmax,it90max] = max(max_s_spec(Nispec:Nfspec));
        Nspec90 = Nispec + it90max;
        Nispec = Nispec+1;
        Nfspec = Nfspec-1;
    end
    t90 = (Nspec90*tempoTotal)/length(max_s_spec);
%Terminava aqui
%end

% N_spec_win = (N_win*length(max_s_spec))/length(x1);
% media_s_spec = zeros(1, length(s_spec(:,Nspec90))); 
% for j=1:length(s_spec(:,Nspec90))
%     media_s_spec(j) = mean(s_spec(j,(Nspec90 - N_spec_win/2):(Nspec90 + N_spec_win/2)));
% end
% media_s_spec = media_s_spec/max(media_s_spec);
% f_spec = f_spec/1000;