function [y_inf, y_sup, t90f] = image_fit (t, t90, threashold, Cmat, tau, v, sentido, plotar, outlier_option)
% Duas buscas. De cima para baixo e de baixo para cima.

% Entradas:
% t - vetor de tempo
% t90 - instante da passagem do carro em frente ao array 
% threashold - limiar para eliminar dados irrelevantes
% Cmat - matriz contendo a funcao da qual se deseja obter os picos
% tau - vetor dos atrasos possiveis
% plotar - flag para plotar resultados
% outlier_option - `exclude` para eliminar dados discrepantes ou `replace`
%                   para substituir pelos dados da primeira curva ajustada

    % Setando valores default
    if isempty(threashold)
        threashold = 0.08;
    end
    
    if isempty(outlier_option)
        outlier_option = 'exclude';
    end
    
    %tau = tau/1000;
    
    % Calculando janela de interesse (ignora dados distantes de t90)
    t_min = t90 - 1;
    t_max = t90 + 1;
    delta_t = t(end) / (length(t) - 1); % Tempo entre amostras de tempo
    indice_min = floor(t_min / delta_t) + 1; % Indice correspondente ao t_min  
    indice_max = floor(t_max / delta_t) + 1; % Indice correspondente ao t_max
    if indice_min < 1
        indice_min = 1;
    elseif indice_max > length(t)
        indice_max = length(t);
    end
    window = indice_min : indice_max; % Intervalo da janela
    tw = t(window);

    % Aplicando threashold e janela
    Cmat_window = Cmat;
    Cmat_window(Cmat_window < threashold) = 0; % Aplica threashold
    
    close
    image(100*Cmat_window)
    set(gca,'YDir','reverse');
    pause
    
    Cmat_window = (Cmat_window(:,window));   % Extrai janela de tempo de interesse de Cmat
    
    
    
    curva_sup = [0];
    curva_inf = [0];
    
    % Busca pelo picos em cada coluna
    for coluna = 1:length(window)
        [peaks, locs, largura] = findpeaks(Cmat_window(:,coluna)); % Encontra valores de Cmat_window == 1
               
        for i_peak = length(peaks) : -1 : 1  % Iteração iniciando do maior indice para avaliar a curva superior    
            if largura(i_peak) > 1 & length(curva_sup) < coluna % Pelo menos 3 pontos "agrupados", para evitar ruído. E coluna atual não preenchida
                curva_sup(coluna) = locs(i_peak); % Guarda ponto encontrado 
            end
        end

        if length(curva_sup) < coluna % Quando não identificar nenhum pico válido em um instante de tempo, copia o anterior
            curva_sup(coluna) = curva_sup(coluna - 1);
        end

        for i_peak = 1 : length(peaks) % Iteração iniciando do menor indice para avaliar a curva inferior
            if largura(i_peak) > 0  & length(curva_inf) < coluna % Pelo menos 2 pontos "agrupados", para evitar ruído. E coluna atual não preenchida
                curva_inf(coluna) = locs(i_peak); 
            end
        end

        if length(curva_inf) < coluna % Quando não identificar nenhum pico válido em um instante de tempo, copia o anterior
            curva_inf(coluna) = curva_inf(coluna - 1);
        end

    end
    
    for coluna = 2:(length(curva_inf)-1) % Inicia na segunda coluna para sempre ter o ponto anterior

        % Condição abaixo verifica se o ponto da curva inferior está
        % afastado a um intervalo de tolerância do ponto da curva
        % superior.
       
        tolerancia = abs(max(curva_sup) - min(curva_sup))*0.08;
        if (abs(curva_sup(coluna) - curva_inf(coluna)) < tolerancia)
            curva_inf(coluna) = curva_inf(coluna-1);
        end       
        
    end
    
    % Encontrando ajustes das curvas
    curva_sup = (((curva_sup./max(curva_sup)).*(2*max(tau))) + min(tau));
    curva_inf = (((curva_inf./max(curva_inf)).*(2*max(tau))) + min(tau));
    
    % Verifica sentido de movimento do veículo
    if sentido == '0° -> 180°'
        sinal = -1; 
    else
        sinal = 1;
    end
    
    x = t(window);
    dist_source_mic = 3.6;
    dist_source_source = 2.5;
    dist_mic_mic = 0.2;
    heigth_mic = 1.26;
        
    options = fitoptions('Method', 'NonlinearLeastSquares', 'Startpoint', [v, t90, dist_source_mic],  'Lower',[0, 0, 0], 'Upper', [100, t(end), 10]);
    
    % Fitting da curva superior
    f_sup = fittype([ num2str(sinal),'*(sqrt( (((x - c) .* v/3.6) +', num2str(dist_mic_mic), ').^2 + (a)^2 + (', num2str(heigth_mic),')^2) - sqrt( ((x - c) .* v/3.6 ).^2 + (a)^2 + (', num2str(heigth_mic),')^2 )) / 343'],...
          'dependent', {'y_sup'}, 'independent', {'x'}, 'coefficients', {'v','c','a'});
    fit_sup = fit(t(window).', curva_sup.', f_sup, options);
    
    % Atualizando parâmetros iniciais para o ajuste de curva com os resultados da
    % curva superior e adicionando o atraso entre as fontes
    options = fitoptions('Method', 'NonlinearLeastSquares', 'Startpoint', [fit_sup.c, sinal*dist_source_source],  'Lower', [0, -10], 'Upper', [t(end), 10]);
    
    % Fitting da curva inferior
    f_inf = fittype([num2str(sinal),'*(sqrt( (((x - c) .*',num2str(fit_sup.v),'/3.6) + ', num2str(dist_mic_mic), ' - shift ).^2 + (', num2str(fit_sup.a),')^2 + (', num2str(heigth_mic),')^2) - sqrt( ((x - c) .*', num2str(fit_sup.v),'/3.6 - shift).^2 + (', num2str(fit_sup.a), ')^2 + (', num2str(heigth_mic),')^2 )) / 343'],...
          'dependent', {'y_inf'}, 'independent', {'x'}, 'coefficients', {'c', 'shift'});
    fit_inf = fit(t(window)', curva_inf', f_inf, options);
    
    % Salva curvas para plot final
    y_sup_inicial = fit_sup(t);
    y_inf_inicial = fit_inf(t);
    curva_sup_inicial = curva_sup;
    curva_inf_inicial = curva_inf;
    
    % Melhorando os pontos utilizados no ajuste das curvas
    
    tolerancia = 2e-4;
    %tolerancia = 1.5*std(curva_sup_inicial);
    y_sup_win = y_sup_inicial(window);
    y_inf_win = y_inf_inicial(window);
    
    % Guarda indice dos pontos onde o erro entre a curva ajustada e os
    % dados excede a tolerância
    ind_erro_sup = abs(y_sup_win' - curva_sup) > tolerancia;
    ind_erro_inf = abs(y_inf_win' - curva_inf) > tolerancia; 
    
    switch outlier_option
    
        case 'replace'
            
            % Substitui os dados pelo valor da curva nos pontos outliers
            curva_sup(ind_erro_sup) = y_sup_win(ind_erro_sup);
            curva_inf(ind_erro_inf) = y_inf_win(ind_erro_inf);
            % Inicializa outliers com zeros (logicos)
            outliers_sup = false(1, length(curva_sup));
            outliers_inf = false(1, length(curva_inf));
            
        case 'exclude'
        
            outliers_sup = excludedata(tw, curva_sup, 'indices', ind_erro_sup);
            outliers_inf = excludedata(tw , curva_inf, 'indices', ind_erro_inf);
    
        otherwise
            display('ERRO. Opção não disponível. Escolha entre replace ou exclude.');
            return;
            
    end
    
    % Refazendo os ajustes com os novos pontos
    options = fitoptions('Method', 'NonlinearLeastSquares', 'Startpoint',...
                        [v, t90, dist_source_mic],...
                        'Lower',[0, 0, 0], 'Upper', [100, t(end), 10],...
                        'Exclude', outliers_sup);
   
    % Fitting da curva superior
    f_sup = fittype([num2str(sinal),'*(sqrt( (((x - c) .* v/3.6) +', num2str(dist_mic_mic), ').^2 + (a)^2 + (', num2str(heigth_mic),')^2) - sqrt( ((x - c) .* v/3.6).^2 + (a)^2 + (', num2str(heigth_mic),')^2 )) / 343'],...
          'dependent', {'y_sup'}, 'independent', {'x'}, 'coefficients', {'v','c','a'});
    [fit_sup,gof] = fit(t(window)', curva_sup', f_sup, options);
    t90f = fit_sup.c;
    %writetable(struct2table(gof), [filename, '.txt'])
    %fid = fopen([filename, '.txt'], 'wt');
    %fprintf(fid, gof);
    
    % Atualizando parâmetros iniciais para o ajuste de curva com os resultados da
    % curva superior e adicionando o atraso entre as fontes
    options = fitoptions('Method', 'NonlinearLeastSquares', 'Startpoint', [fit_sup.c, sinal*dist_source_source],...  
                         'Lower', [0, -10], 'Upper', [t(end), 10], 'Exclude', outliers_inf);
    
    % Fitting da curva inferior
    f_inf = fittype([num2str(sinal),'*(sqrt( (((x - c) .*',num2str(fit_sup.v),'/3.6) + ', num2str(dist_mic_mic), ' - shift).^2 + (', num2str(fit_sup.a),')^2 + (', num2str(heigth_mic),')^2) - sqrt( ((x - c) .*', num2str(fit_sup.v),'/3.6 - shift).^2 + (', num2str(fit_sup.a), ')^2 + (', num2str(heigth_mic),')^2 )) / 343'],...
          'dependent', {'y_inf'}, 'independent', {'x'}, 'coefficients', {'c', 'shift'});
    fit_inf = fit(t(window)', curva_inf', f_inf, options);
        
    y_sup = fit_sup(t);
    y_inf = fit_inf(t);
    
    % Alternativa para o ajuste da curva inferior: Aplica atraso na curva 
    % superior
%    fit_shift = fit_sup;
%    fit_shift.shift = dist_source_source;
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% PLOTANDO AJUSTES DE CURVA %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if plotar == true
        close all;
        figure
        img_name = strcat(num2str(v), 'dados&estimativa_inf');

        % Plota os dados
        p_dinf = plot(tw, 1000*curva_inf_inicial, '*');hold on;
        p_dsup = plot(tw, 1000*curva_sup_inicial, '*'); hold on;
        p_dinf.Color = RGB('light blue'); p_dinf.MarkerSize = 4;
        p_dsup.Color = RGB('light orange'); p_dsup.MarkerSize = 4;
        axis ([0 t(end) -0.8 0.8]);
        xlabel('Tempo (s)'), ylabel('Atraso \tau (s)');
        h = [p_dsup p_dinf];
        legend(h, 'Dados Fonte 1', 'Dados Fonte 2',...
                'Location', 'northoutside', 'Orientation', 'horizontal')
%         img_name = strcat(filename, '_DATA_ANTES');
%         print(gcf,img_name,'-dpng','-r300'); % Salvando imagem em arquivo png
%         print(gcf,img_name,'-depsc','-r300'); % Salvando imagem em arquivo eps

        % Plota as curvas
        p_cinf = plot(t, 1000*y_inf_inicial); hold on;
        p_csup = plot(t, 1000*y_sup_inicial);
        p_cinf.Color = RGB('dark blue'); pcinf.LineWidth = 1;
        p_csup.Color = RGB('orange'); p_csup.LineWidth = 1;
        axis ([0 t(end) -0.8 0.8]);
        xlabel('Tempo t (s)'), ylabel('Atraso \tau (s)');
        h = [p_csup p_dsup p_cinf p_dinf];
        legend(h, 'Curva Fonte 1', 'Curva Fonte 2', 'Dados Fonte 1', 'Dados Fonte 2',...
                'Location', 'northoutside', 'Orientation', 'horizontal')
%         img_name = strcat(filename, '_CURVA_ANTES');
%         print(gcf,img_name,'-dpng','-r300'); % Salvando imagem em arquivo png
%         print(gcf,img_name,'-depsc','-r300'); % Salvando imagem em arquivo eps

        % Depois do PÓS-PROCESSAMENTO

        % Plota os dados
        close all;

        if outlier_option == 'replace'
            p_dinf = plot(tw, 1000*curva_inf, '*');hold on;
            p_dsup = plot(tw, 1000*curva_sup, '*'); hold on;
        elseif outlier_option == 'exclude'
            p_dinf = plot(tw(~outliers_inf), 1000*curva_inf(~outliers_inf), '*');hold on;
            p_dsup = plot(tw(~outliers_sup), 1000*curva_sup(~outliers_sup), '*'); hold on;
        end
        p_dinf.Color = RGB('light blue'); p_dinf.MarkerSize = 4;
        p_dsup.Color = RGB('light orange'); p_dsup.MarkerSize = 4;
        axis ([0 t(end) -0.8 0.8]);
        xlabel('Tempo t (s)'), ylabel('Atraso \tau (s)');
        h = [p_dsup p_dinf];
        legend(h, 'Dados Fonte 1', 'Dados Fonte 2',...
                'Location', 'northoutside', 'Orientation', 'horizontal')

        % Plota as curvas
        p_cinf = plot(t, 1000*y_inf); hold on;
        p_csup = plot(t, 1000*y_sup);
        p_cinf.Color = RGB('dark blue'); pcinf.LineWidth = 1;
        p_csup.Color = RGB('orange'); p_csup.LineWidth = 1;
        axis ([0 t(end) -0.8 0.8]);
        xlabel('Tempo t (s)'), ylabel('Atraso \tau (s)');
        h = [p_csup p_dsup p_cinf p_dinf];
        legend(h, 'Curva Fonte 1', 'Curva Fonte 2', 'Dados Fonte 1', 'Dados Fonte 2',...
                'Location', 'northoutside', 'Orientation', 'horizontal')
        
    end
    
end        