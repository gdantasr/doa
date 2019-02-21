function [y_inf, y_sup, fit_inf, fit_sup] = edge_detect (t, t90, tau, Cmat, v, sentido, plotar, outlier_option, filename, dist)
%
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
    
    % Normaliza cada coluna da matriz de correlação
    Cmat = Cmat / max(abs(Cmat(:)));
            
    % Calculando janela de interesse (ignora dados distantes de t90)
    t_min = t90 - 1;
    t_max = t90 + 1;
    
    [~, indice_min] = min(abs(t_min - t));
    [~, indice_max] = min(abs(t_max - t));

    window = indice_min : indice_max; % Intervalo da janela
    tw = t(window);
    Cwindow = Cmat(:, window);
    
    threshold = 0.15*max(Cwindow(:));
    
    G = mat2gray(Cwindow, [threshold max(Cwindow(:))]);
%     image(100*G)

    %1
    se = strel('disk', 1);
    BW = imopen(G, se);
%     image(100*BW)
        
    %8
    se = strel('disk', 8);
    BW = imclose(BW, se);
%     image(100*BW)
    
    %5
    se = strel('disk', 5);
    BW = imopen(BW, se);
%     image(100*BW)
               
    curva_sup = [0];
    curva_inf = [0];
    cpeaks = zeros(size(BW));
    for coluna = 1:size(BW,2),
        [peaks, locs, largura] = findpeaks(BW(:,coluna));
        tabPicos = table(locs, peaks, largura, 'VariableNames', {'Posicao' ; 'Amplitude' ; 'Largura'} );
        tabPicos = sortrows(tabPicos, 'Amplitude');
        
        if length(peaks) == 0 % Nenhum pico encontrado
             [~, mainPeaks] = max(BW(:,coluna));
             curva_sup(coluna) = mainPeaks(1);
             curva_inf(coluna) = mainPeaks(1);
        elseif length(peaks) == 1 % Apenas 1 pico encontrado
             mainPeaks = tabPicos.Posicao(1);
             curva_sup(coluna) = mainPeaks(1);
             curva_inf(coluna) = mainPeaks(1);
        elseif (length(peaks) > 2 || length(peaks) == 2)
            mainPeaks = sort(tabPicos.Posicao(1:2), 'descend'); % Posição dos 2 picos principais em ordem decrescente da posição
            curva_sup(coluna) = mainPeaks(1);
            curva_inf(coluna) = mainPeaks(2);
        end
        cpeaks (mainPeaks, coluna) = 1;

    end
    
    % Encontrando ajustes das curvas
    curva_sup = (((curva_sup./max(curva_sup)).*(2*max(tau))) + min(tau));
    curva_inf = (((curva_inf./max(curva_inf)).*(2*max(tau))) + min(tau));
    
    % Verifica sentido de movimento do veículo
    if strcmp (sentido, '0° -> 180°')
        sinal = -1; 
    else
        sinal = 1;
    end
    
    x = t(window);
    vs = 340;
        
    options = fitoptions(   'Method',       'NonlinearLeastSquares', ... 
                            'Startpoint',   [t90, dist.source_mic, v], ... 
                            'Lower',        [0, 0, 0], ...
                            'Upper',        [t(end), 10, 100]);   
  
    % Fitting da curva superior
    f_sup = fittype(@(a,b,v,x)passby( a, b, v, x, dist.mic_mic, dist.heigth, sinal), ...
           'independent', {'x'}, 'coefficients', {'a','b','v'});
    fit_sup = fit(t(window)', curva_sup', f_sup, options);
    
    % Atualizando parâmetros iniciais para o ajuste de curva com os resultados da
    % curva superior e adicionando o atraso entre as fontes
    options = fitoptions(   'Method',       'NonlinearLeastSquares', ...
                            'Startpoint',   [fit_sup.a], ... 
                            'Lower',        [0], ... 
                            'Upper', [t(end)]);
    
    % Fitting da curva inferior
    f_inf = fittype( @(a,x)passby(a, fit_sup.b, fit_sup.v, x, dist.mic_mic, dist.heigth, sinal), ...
              'independent', {'x'}, 'coefficients', {'a'} );
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
    options = fitoptions('Method', 'NonlinearLeastSquares', ...
                        'Startpoint', [t90, dist.source_mic, v],...
                        'Lower',[0, 0, 0], 'Upper', [t(end), 10, 100],...
                        'Exclude', outliers_sup);
   
    % Fitting da curva superior
    %[fit_sup, gof] = fit(...);
    fit_sup = fit(t(window)', curva_sup', f_sup, options);
    t90f = fit_sup.a;
    
    %writetable(struct2table(gof), [filename, '.txt'])
    %fid = fopen([filename, '.txt'], 'wt');
    %fprintf(fid, gof);
    
    % Atualizando parâmetros iniciais para o ajuste de curva com os resultados da
    % curva superior e adicionando o atraso entre as fontes
    options = fitoptions(   'Method',       'NonlinearLeastSquares', ...
                            'Startpoint',   [fit_sup.a], ...  
                            'Lower',        [0], ...
                            'Upper',        [t(end)], ...
                            'Exclude', outliers_inf);  % NAO TAVA USANDO ESSE EXCLUDE. POR QUE?
                         
    % Fitting da curva inferior
%     f_inf = fittype([num2str(sinal),'*(sqrt( (((x - c) .*',num2str(fit_sup.v),'/3.6) + ', num2str(dist.mic_mic), ').^2 + (', num2str(fit_sup.a),')^2 + (', num2str(dist.heigth),')^2) - sqrt( ((x - c) .*', num2str(fit_sup.v),'/3.6).^2 + (', num2str(fit_sup.a), ')^2 + (', num2str(dist.heigth),')^2 )) / 343'],...
%           'dependent', {'y_inf'}, 'independent', {'x'}, 'coefficients', {'c'});
    f_inf = fittype( @(a,x)passby(a, fit_sup.b, fit_sup.v, x, dist.mic_mic, dist.heigth, sinal), ...
                     'independent', {'x'}, 'coefficients', {'a'} );
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
        xlabel('Tempo (s)'), ylabel('Atraso \tau (ms)');
        h = [p_dsup p_dinf];
        legend(h, 'Dados Fonte 1', 'Dados Fonte 2',...
                'Location', 'northoutside', 'Orientation', 'horizontal')
        img_name = strcat(filename, '_DATA_ANTES');
        print(gcf,img_name,'-dpng','-r300'); % Salvando imagem em arquivo png
        print(gcf,img_name,'-depsc','-r300'); % Salvando imagem em arquivo eps

        % Plota as curvas
        p_cinf = plot(t, 1000*y_inf_inicial); hold on;
        p_csup = plot(t, 1000*y_sup_inicial);
        p_cinf.Color = RGB('dark blue'); pcinf.LineWidth = 1;
        p_csup.Color = RGB('orange'); p_csup.LineWidth = 1;
        axis ([0 t(end) -0.8 0.8]);
        xlabel('Tempo t (s)'), ylabel('Atraso \tau (ms)');
        h = [p_csup p_dsup p_cinf p_dinf];
        legend(h, 'Curva Fonte 1', 'Curva Fonte 2', 'Dados Fonte 1', 'Dados Fonte 2',...
                'Location', 'northoutside', 'Orientation', 'horizontal')
        img_name = strcat(filename, '_CURVA_ANTES');
        print(gcf,img_name,'-dpng','-r300'); % Salvando imagem em arquivo png
        print(gcf,img_name,'-depsc','-r300'); % Salvando imagem em arquivo eps

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
