clear all; clc; close all;

%% Definição dos parâmetros
vs = 340; % Velocidade de propagação do som
d= 0.2;
dist_source_mic = 3.6;
dist_source_source = 2.5;
dist_mic_mic = 0.2;
heigth_mic = 1.26;
mic1=1;   % mic1 e mic2 são os microfones
mic2=2;   % usados nos algorítimos de estimação do DOA
Fresample = 32000; % Alterar caso queira reamostrar o sinal
N = 256;
plot_doa = true;  % Flag para plotar a curva da DOA
plot_fit = true; % Flag para plotar curvas ajustadas (considerando as duas fontes)
mics=[-0.2 0; 0 0; 0.2 0; 0 -0.2]; %microphone array geometry

pathAudio = strcat('/home/netware/users/gabidantas/Documents/Mestrado/TEPS/Audios/'); % Caminho para a pasta com os áudios
pathAudioCut = strcat('/home/netware/users/gabidantas/Documents/Mestrado/TEPS/BancoAudios/'); % Caminho para a pasta com os áudios
for indiceTeste = 1:10
    
    %% Lendo arquivo de texto
    fileName = strcat(pathAudioCut, 'Teste', num2str(indiceTeste), '.txt'); % Nome do arquivo txt
    fileID = fopen(fileName,'r');

    fgetl(fileID); % Pula as duas
    fgetl(fileID); % primeiras linhas

    velocidade = fgetl(fileID); % Extrai linha que contem velocidade
    sentido = fgetl(fileID);    % Extrai linha que contem sentido
    qtdeCarros = fgetl(fileID); % Extrai linha que contem qtde de carros

    % Ignorando as próximas linhas 
    fgetl(fileID); fgetl(fileID); fgetl(fileID); fgetl(fileID);

    instantes = fgetl(fileID);

    fclose(fileID);

    % Extraindo dados das linhas do arquivo de texto
    velocidade = regexp(velocidade, '\d*', 'match');
    velocidade = str2num(char(velocidade));
    sentido = sentido(10:19);
    qtdeCarros = regexp(qtdeCarros, '\d*', 'match');
    qtdeCarros = str2num(char(qtdeCarros));
    
    % Lê audios originais, que contém trechos "de silêncio"
    fileName = strcat(pathAudio,'Teste', num2str(indiceTeste), '/1 Faixa de Audio.wav'); % Lendo apenas mic 1
    [noise, Fs] = audioread(fileName);
    noise = noise(1:Fs); % Trecho do sinal onde há só ruído

    %% Lendo arquivo de audio
    for indiceCarro = 1:qtdeCarros
        Yxt = [];
        for indiceMic = 1:4
            fileName = strcat(pathAudioCut,'teste', num2str(indiceTeste),'_carro', num2str(indiceCarro) ,'_v', num2str(velocidade),'_mic',num2str(indiceMic),'.wav');
            [Yx,Fs] = audioread(fileName);
            Yxt(1:length(Yx),indiceMic)=Yx/max(Yx);
        end

        %% Calculando DOA
        
        % Pré-processamento do sinal: Filtragem e reamostragem
        Fm = 200; % Frequencia de corte inferior
        Fc = 4000; % Frequencia de corte do filtro
        [b,a] = butter(5, [Fm/(Fs/2) Fc/(Fs/2)], 'bandpass'); 
        Yxt_filt = filter(b, a, Yxt); % Filtra
        Yxt_resample = resample(Yxt_filt, Fresample, Fs); % Reamostra
        
        % Calcula instante da passagem do veículo t90
        n90 = energy_peak(Yxt_resample(:,1), Fresample);
        t90 = n90/Fresample;

        % Método GCC-PHAT para estimar a DOA       
        [phi, tdd, t, tau, Cmat] = doa_gcc(Yxt_resample(:,mic1), Yxt_resample(:,mic2), 0.2, N, Fresample);
        
%         % Normaliza cada coluna da matriz de correlação
%         Cmatn = zeros(size(Cmat));
%         for indice = 1:size(Cmat,2)
%            Cmatn(:, indice) = Cmat(:, indice) ./ max(Cmat(:,indice));
%         end
%         Cmat = Cmatn;
        
        % Threshold da matriz de correlação
        th = 0.15*max(max(Cmat));
        
        % Pós-processamento
        filename = strcat('TDD_teste', num2str(indiceTeste),'_carro', num2str(indiceCarro) ,'_v', num2str(velocidade), '_f', num2str(Fresample),'_N', num2str(N), '_th', num2str(th));
        
        [tdd_inf, tdd_sup, fit_inf, fit_sup] = edge_detect (t, t90, th, tau, Cmat, velocidade, sentido, plot_fit, 'exclude', filename);
        %[tdd_inf, tdd_sup, t90f] = image_fit (t, t90, th, Cmat, tau, velocidade, sentido, plot_fit, 'exclude');
        
%         % Salva parâmetros do ajuste para usar na curva teórica
%         v = fit_sup.v;
         t90_sup = fit_sup.c;
         t90_inf = fit_inf.c;
%         dist_source_mic = fit_sup.a;
        
        phi_inf = 180/pi*real(acos(vs/d*tdd_inf)); phi_sup = 180/pi*real(acos(vs/d*tdd_sup));
        
        %titulo =  strcat('TDD - Fs=', num2str(Fresample),'Hz, N=', num2str(N),' Lov=25%, OV=4, \alpha=0.8');
        titulo =  strcat('Função GCC-PHAT');
        plot_tdd(t, tau, [tdd_sup tdd_inf], Cmat, false, filename, titulo);
        %pause
        
        % Teórico
%        [phi_theo_1, phi_theo_2, tdd_theo, tdd_theo_2] = phi_theo_multisource (velocidade, t90_inf, t90_sup, sentido, dist_source_mic, dist_source_source, heigth_mic, T);
%        phi_resample = resample(phi, length(phi_theo_1), length(phi));

        % Beamformer
        [y_sup, BP] = adapt_beam(Yxt_resample, mics, phi_sup, 0.01, 4, 1000, Fresample);
        [y_inf, BP] = adapt_beam(Yxt_resample, mics, phi_inf, 0.01, 4, 1000, Fresample);
        
 %       audiowrite(['./beamformer/beam_sup_t', num2str(indiceTeste), '_c', num2str(indiceCarro) ,'_v', num2str(velocidade), '_f', num2str(Fresample),'_N', num2str(N), '_th', num2str(th), '.wav'], y_sup, Fresample);
 %       audiowrite(['./beamformer/beam_inf_t', num2str(indiceTeste), '_c', num2str(indiceCarro) ,'_v', num2str(velocidade), '_f', num2str(Fresample),'_N', num2str(N), '_th', num2str(th), '.wav'], y_inf, Fresample);
        
        % Interpolando
 %       audiowrite(['./beamformer/rs_beam_sup_t', num2str(indiceTeste), '_c', num2str(indiceCarro) ,'_v', num2str(velocidade), '_f', num2str(Fresample),'_N', num2str(N), '_th', num2str(th), '.wav'], resample(y_sup, Fs, Fresample), Fs);
 %       audiowrite(['./beamformer/beam_inf_t', num2str(indiceTeste), '_c', num2str(indiceCarro) ,'_v', num2str(velocidade), '_f', num2str(Fresample),'_N', num2str(N), '_th', num2str(th), '.wav'], resample(y_inf, Fs, Fresample), Fs);
        
        % Curva Superior
        figure
        subplot(3,1,1)
        spectrogram(Yxt_resample(:,1) ,hamming(N), N/2, N, Fresample,'yaxis');
        colormap('jet');
        title(['Espectro da Entrada do Beamforming. Fs = ', num2str(Fresample),'Hz']);
        xlabel('Tempo(s)')
        ylabel('Frequência(Hz)')
        
        subplot(3,1,2)
        [s_sup, w_sup, t_sup] = spectrogram(y_sup ,hamming(N), N/2, N, Fresample,'yaxis');
        
        i_spec = 1;
        while t_sup(i_spec) < t90_sup
            i_spec = i_spec + 1;
        end
        plot(w_sup, 10*log10(abs(s_sup(:,i_spec))))
        ylim([-40 10])
        title(['Espectro da Saída do Beamforming (Roda Dianteira). Fs = ', num2str(Fresample),'Hz']);
        
        subplot(3,1,3)
        [s_inf, w_inf, t_inf] = spectrogram(y_inf ,hamming(N), N/2, N, Fresample,'yaxis');
        
        i_spec = 1;
        while t_sup(i_spec) < t90_sup
            i_spec = i_spec + 1;
        end
        plot(w_inf, 10*log10(abs(s_inf(:,i_spec))))
        ylim([-40 10])
        title(['Espectro da Saída do Beamforming (Roda Traseira). Fs = ', num2str(Fresample),'Hz']);
        
        fig = gcf;
        fig.InvertHardcopy = 'off';
        
        filename = ['./beamformer/beam_spec_t', num2str(indiceTeste), '_c', num2str(indiceCarro) ,'_v', num2str(velocidade), '_f', num2str(Fresample),'_N', num2str(N), '_th', num2str(th)];
%        print(gcf, filename,'-dpng','-r300'); % Salvando imagem em arquivo png

%         filename = strcat('beam_spec_fs',num2str(Fresample), '_f', num2str(fbeam),'.png');
%         print(gcf, filename,'-dpng','-r300'); % Salvando imagem em arquivo png
        
%         % Curva inferior
%         figure
%         subplot(2,1,1)
%         spectrogram(yx_resample(:,1) ,hamming(N), N/2, N, Fresample,'yaxis');
%         colormap('jet');
%         title(['Beamforming Input Signal Spectrogram. Fs = ', num2str(Fresample),'Hz']);
%         subplot(2,1,2)
%         spectrogram(y_inf ,hamming(N), N/2, N, Fresample,'yaxis');
%         colormap('jet');
%         title(['Beamforming Output Signal Spectrogram. Fs = ', num2str(Fresample),'Hz']);
%         fig = gcf;
%         fig.InvertHardcopy = 'off';        
% %         filename = strcat('beam_spec_fs',num2str(Fresample), '_f', num2str(fbeam),'.png');
% %         print(gcf, filename,'-dpng','-r300'); % Salvando imagem em arquivo png
        
    end
end

close all;