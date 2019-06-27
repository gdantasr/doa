function [y_inf, y_sup, fit_inf, fit_sup] = fit_curve2mat (t, t90, tau, Cmat, v, sentido, plotar, outlier_option, filename, dist, th)
%
% Duas buscas. De cima para baixo e de baixo para cima.

% Entradas:
% t - vetor de tempo
% t90 - instante da passagem do carro em frente ao array
% tau - vetor dos atrasos possiveis
% Cmat - matriz contendo a funcao da qual se deseja obter os picos
% plotar - flag para plotar resultados
% filename - caso queira salvar os plots
% outlier_option - `exclude` para eliminar dados discrepantes ou `replace`
%                   para substituir pelos dados da primeira curva ajustada
% dist - estrutura contento as distâncias do problema
    
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
    threshold = th*max(Cwindow(:)); 
    
    %% "Image" processing steps
    
    % Structuring elements size
    t_s = mean(diff(t));        % Time sampling period
    xWidth = dist.wheelbase * 3.6 / (v * t_s);
    yWidth = 4;
    se_dilate = strel('rectangle', [round(xWidth)+1 yWidth]);
    se_erode = strel('rectangle', [round(xWidth)-1 yWidth]);
    
    % Scale CrossCorr matrix values to grayscale
    img_gray = mat2gray(Cwindow, [threshold max(Cwindow(:))]);
    % Dilatation
    img_dilated = imdilate(img_gray, se_dilate);
    % Erosion
    img_eroded = imerode(img_dilated, se_erode);
    % Grayscale to Binary
    img_binary = im2bw(img_eroded, th*max(img_eroded(:)));
    % "Remove" morphologic operation 
    img = bwmorph(img_binary,'remove');

    %% Get data from image
    
    % Find mean curve from image
    mean_curve = zeros(size(window));
    for col = 1:size(img, 2)
        mean_curve(col) = round(mean(find(img(:, col))));
    end
    % Interpolate to remove NaN data
    xdata = (1:length(mean_curve))';
    mean_curve = interp1(xdata(~isnan(mean_curve)), mean_curve(~isnan(mean_curve)), xdata, 'linear', 'extrap');
    over = mean_curve > size(img,1);
    mean_curve(over) = size(img,1);
    over = mean_curve < 1;
    mean_curve(over) = 1;
    % Median filter to remove noisy peaks
    mean_curve = round(medfilt1(mean_curve,3));
        
    %% Prepare data for curve fitting
    
    % Separate up and bottom curves by comparing with mean curve
    curva_sup = zeros(size(mean_curve'));
    curva_inf = zeros(size(mean_curve'));
    for col = 1:size(img,2)                 % Loop through image columns
        non_zero_lines = find(img(:, col)); % Find nonzero line indexes
        is_sup = non_zero_lines > mean_curve(col);
        % Top curve
        if ~isempty(find(is_sup, 1))         % If there are any nonzero lines...
            curva_sup(col) = tau(min(non_zero_lines(is_sup)));
        else
            curva_sup(col) = tau(mean_curve(col));
        end
        % Bottom curve
        if ~isempty(find(~is_sup, 1))         % If there are any nonzero lines...
            curva_inf(col) = tau(max(non_zero_lines(~is_sup))); 
        else
            curva_inf(col) = tau(mean_curve(col));
        end
    end

    % Check pass-by direction
    if strcmp (sentido, '0° -> 180°')
        sinal = -1; 
    else
        sinal = 1;
    end
    
    % Define options for first fitting step
    options = fitoptions(   'Method',       'NonlinearLeastSquares', ... 
                            'Startpoint',   [t90, dist.source_mic, v], ... 
                            'Lower',        [0, 0, 0], ...
                            'Upper',        [t(end), 10, 100]);   
  
    % Upper curve fitting
    f_sup = fittype(@(a,b,v,x)passby( a, b, v, x, dist.mic_mic, dist.heigth, sinal), ...
           'independent', {'x'}, 'coefficients', {'a','b','v'});
    
    fit_sup = fit(t(window)', curva_sup', f_sup, options);
    
    % Update startpoints with upper curve fit results
    options = fitoptions(   'Method',       'NonlinearLeastSquares', ...
                            'Startpoint',   [fit_sup.a], ... 
                            'Lower',        [0], ... 
                            'Upper', [t(end)]);
    
    % Bottom curve fitting
    f_inf = fittype( @(a,x)passby(a, fit_sup.b, fit_sup.v, x, dist.mic_mic, dist.heigth, sinal), ...
              'independent', {'x'}, 'coefficients', {'a'} );
    fit_inf = fit(t(window)', curva_inf', f_inf, options);
    
    % Save curves for final plot
    y_sup_inicial = fit_sup(t);
    y_inf_inicial = fit_inf(t);
    curva_sup_inicial = curva_sup;
    curva_inf_inicial = curva_inf;
    
    % Improve fitting    
    tolerancia = 2e-4;
    %tolerancia = 1.5*std(curva_sup_inicial);
    y_sup_win = y_sup_inicial(window);
    y_inf_win = y_inf_inicial(window);
    
    % Find outliers
    ind_erro_sup = abs(y_sup_win' - curva_sup) > tolerancia;
    ind_erro_inf = abs(y_inf_win' - curva_inf) > tolerancia; 
    
    switch outlier_option
    
        case 'replace'
            
            % Substitute outliers with initial fit curve points
            curva_sup(ind_erro_sup) = y_sup_win(ind_erro_sup);
            curva_inf(ind_erro_inf) = y_inf_win(ind_erro_inf);
            % Outliers vector initialization
            outliers_sup = false(1, length(curva_sup));
            outliers_inf = false(1, length(curva_inf));
            
        case 'exclude'
        
            outliers_sup = excludedata(tw, curva_sup, 'indices', ind_erro_sup);
            outliers_inf = excludedata(tw , curva_inf, 'indices', ind_erro_inf);
    
        otherwise
            display('ERRO. Opção não disponível. Escolha entre replace ou exclude.');
            return;
    end
    
    % Repeat fitting step with new data points (no outliers)
    options = fitoptions('Method', 'NonlinearLeastSquares', ...
                        'Startpoint', [t90, dist.source_mic, v],...
                        'Lower',[0, 0, 0], 'Upper', [t(end), 10, 100],...
                        'Exclude', outliers_sup);
   
    % Upper curve fitting
    %[fit_sup, gof] = fit(...);
    fit_sup = fit(t(window)', curva_sup', f_sup, options);
    t90f = fit_sup.a;
    
    %writetable(struct2table(gof), [filename, '.txt'])
    %fid = fopen([filename, '.txt'], 'wt');
    %fprintf(fid, gof);
    
    % Update startpoints with upper curve fit results
    options = fitoptions(   'Method',       'NonlinearLeastSquares', ...
                            'Startpoint',   [fit_sup.a], ...  
                            'Lower',        [0], ...
                            'Upper',        [t(end)], ...
                            'Exclude', outliers_inf);  % NAO TAVA USANDO ESSE EXCLUDE. POR QUE?
                         
    % Bottom curve fitting
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