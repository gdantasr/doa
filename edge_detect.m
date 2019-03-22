function [y_inf, y_sup, fit_inf, fit_sup] = edge_detect (t, t90, tau, Cmat, v, sentido, plotar, outlier_option, filename, dist, theo)
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
    
    % Cross corr matrix normalization
    Cmat = Cmat/max(abs(Cmat(:)));
            
    % Time window definition (around t90, when the car crosses the array)
    t_min = t90 - 2;
    t_max = t90 + 2;
    [~, indice_min] = min(abs(t_min - t));
    [~, indice_max] = min(abs(t_max - t));
    window = indice_min : indice_max; % Intervalo da janela
    tw = t(window);
    Cwindow = Cmat(:, window);
    threshold = 0.1*max(Cwindow(:));
    
    %% "Image" processing steps
    
    % Scale CrossCorr matrix values to grayscale
    img_gray = mat2gray(Cwindow, [threshold max(Cwindow(:))]);
%     image(100*G)
    
    % Dilatation
    se = strel('disk', 10);
    img_dilated = imdilate(img_gray, se);
%     image(100*BW);
    
    % Erosion
    se = strel('disk', 2);
    img_eroded = imerode(img_dilated, se);
%     image(100*img_eroded);
    
    % Grayscale to Binary
    img_binary = im2bw(img_eroded, 0.5*max(img_eroded(:)));
%     image(100*img_binary);
   
    % "Remove" morphologic operation 
    img = bwmorph(img_binary,'remove');
%     image(100*img);
    
%     [nonZeroRows, nonZeroCols] = find(BW4); % Pontos não nulos. Pontos das "bordas"
%     img2curve = zeros(size(img_gray)); % Inicializa
%     curve = zeros(size(window));
    
%     i = 1;
%     while i<=length(nonZeroCols)
%         j = i;
%         while j+1 < length(nonZeroCols) && nonZeroCols(j+1)==nonZeroCols(j)
%             j = j+1; % Enquanto j marca a mesma coluna, agrupa
%         end
%         mean_point = round(mean(nonZeroRows(i:j))); % Pega o ponto médio entre as linhas não nulas da mesma coluna
%         img2curve(mean_point, nonZeroCols(i)) = 1;
%         curve(nonZeroCols(i)) = mean_point;
%         i = j+1;
%     end
    
    % Find mean curve from image
    mean_curve = zeros(size(window));
    for col = 1:size(img, 2)
        mean_curve(col) = round(mean(find(img(:, col))));
    end
        
    % Interpolate to remove Nan data
    xdata = (1:length(mean_curve))';
    mean_curve = interp1( xdata(~isnan(mean_curve)), mean_curve(~isnan(mean_curve)), xdata);

    % Median filter to remove noisy peaks
    mean_curve = round(medfilt1(mean_curve,3));
    
    % Encontrando ajustes das curvas
%     curva_sup = (((curva_sup./max(curva_sup)).*(2*max(tau))) + min(tau));
%     curva_inf = (((curva_inf./max(curva_inf)).*(2*max(tau))) + min(tau));
%     curve_scaled = zeros(size(curve));
%     for i = 1:length(curve)
%         if curve(i) ~= 0
%             curve_scaled(i) = tau(curve(i));
%         end
%     end
    

    
%     windowSize = 5; 
%     b = (1/windowSize)*ones(1,windowSize);
%     a = 1;
%     filt_curve = filter(b,a,filt_curve1);
    
    % Verifica sentido de movimento do veículo
    if strcmp (sentido, '0° -> 180°')
        sinal = -1; 
    else
        sinal = 1;
    end
        
    x = t(window);
    vs = 340;
           
    imagesc(tw, tau, img);
    colormap(flipud(bone)); %colorbar;
    set(gca,'YDir','normal'); hold on
    plot(tw, mean_curve, 'b'); hold on
    plot(tw, mean_curve, 'r');
    
    % Separate up and bottom curves by comparing with mean curve
    curva_sup = zeros(size(mean_curve'));
    curva_inf = zeros(size(mean_curve'));
    for col = 1:size(img,2)                 % Loop through image columns
        non_zero_lines = find(img(:, col)); % Find nonzero line indexes
        if ~isempty(non_zero_lines)         % If there are any nonzero lines...
            is_sup = non_zero_lines > mean_curve(col);
            curva_sup(col) = tau(min(non_zero_lines(is_sup)));     % Sup and Inf curves keep lines
            curva_inf(col) = tau(max(non_zero_lines(~is_sup)));    % which are closest to mean curve
        else                                % Else, use mean curve for both
            curva_sup(col) = tau(mean_curve(col));
            curva_inf(col) = tau(mean_curve(col));
        end
    end
       
    % Define options for first fitting step
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