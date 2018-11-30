function F = plot_spec( data, fs, N, speed, carID, plotOption )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Plot options

switch plotOption

    case 0
        %% Plot spectrogram from all sensors
        for fileID = 1 : length(fileNames)
            figure;
            set(gcf, 'Position', get(0, 'Screensize'))
            for i = 1 : 11
                subplot(3,4,i)
                spectrogram(data{fileID}(:,i), hamming(N), N/2, 1024, fs,'MinThreshold', -80 ,'yaxis')
            end
            titleHandle = suptitle(['Espectrograma dos 11 sensores. Teste ', fileNames{fileID}]);
            set(titleHandle,'FontSize',18,'FontWeight','bold')
            colormap(parula(fs));
            F = getframe(gcf);

        end


    case 1
        %% Plot spectrogram: constant speed vs accelerating (full time)       
        figure;
        set(gcf, 'Position', get(0, 'Screensize'))
        
        [~,~,~,P] = spectrogram(data{1}(:,1), hamming(N), N/2, 1024, fs,'yaxis'); % Signal PSD
        specMin = floor(10*log10(min(P(:))));   % Min and max PSD values to 
        specMax = ceil(10*log10(max(P(:))));    % limit colormap range
        
        % Constant speed
        subplot(2,1,1)
        spectrogram(data{1}(:,1), hamming(N), N/2, 1024, fs,'yaxis'); % Plots spectrogram
        colormap(jet(fs));  % Choose plot color scheme
        caxis([-80 -40]);   % Choose colomap range (values outside range are clipped to limits)
        title(['Velocidade constante (', speed{1}, ' km/h)']); % Subplot title
                
        % Accelerating
        subplot(2,1,2)
        spectrogram(data{2}(:,1), hamming(N), N/2, 1024, fs,'yaxis')
        colormap(jet(fs));
        caxis([-80 -40]);
        title('Aceleração');
        
        titleHandle = suptitle(['Espectrograma. Velocidade Constante vs Aceleração. Carro ', carID{1},', sensor 1.']); % Image title
        set(titleHandle,'FontSize',16,'FontWeight','bold')
        F = getframe(gcf);
        
        
        case 2
        %% Plot spectrogram: constant speed vs accelerating (cropped)       
        figure;
        set(gcf, 'Position', get(0, 'Screensize'))
        
        [~,~,~,P] = spectrogram(data{1}(:,1), hamming(N), N/2, 1024, fs,'yaxis'); % Signal PSD
        specMin = floor(10*log10(min(P(:))));   % Min and max PSD values to 
        specMax = ceil(10*log10(max(P(:))));    % limit colormap range
        
        % Constant speed
        subplot(2,1,1)
        spectrogram(data{1}(:,1), hamming(N), N/2, 1024, fs,'yaxis'); % Plots spectrogram
        colormap(jet(fs));  % Choose plot color scheme
        caxis([-80 -40]);   % Choose colomap range (values outside range are clipped to limits)
        title(['Velocidade constante (', speed{1}, ' km/h)']); % Subplot title
                
        % Accelerating
        subplot(2,1,2)
        spectrogram(data{2}(:,1), hamming(N), N/2, 1024, fs,'yaxis')
        colormap(jet(fs));
        caxis([-80 -40]);
        title('Aceleração');
        
        titleHandle = suptitle(['Espectrograma. Velocidade Constante vs Aceleração. Carro ', carID{1},', sensor 1.']); % Image title
        set(titleHandle,'FontSize',16,'FontWeight','bold')
        F = getframe(gcf);       


        case 3
        %% Plot spectrogram: constant speed (all speeds) vs accelerating (fulltime)       
        figure;
        set(gcf, 'Position', get(0, 'Screensize'))
        
        [~,~,~,P] = spectrogram(data{1}(:,1), hamming(N), N/2, 1024, fs,'yaxis'); % Signal PSD
        specMin = floor(10*log10(min(P(:))));   % Min and max PSD values to 
        specMax = ceil(10*log10(max(P(:))));    % limit colormap range
        
        % Constant speed 1
        subplot(4,1,1)
        spectrogram(data{1}(:,1), hamming(N), N/2, 1024, fs,'yaxis'); % Plots spectrogram
        colormap(jet(fs));  % Choose plot color scheme
        caxis([-80 -40]);   % Choose colomap range (values outside range are clipped to limits)
        title(['Velocidade constante (', speed{1}, ' km/h)']); % Subplot title
        
        % Constant speed 2
        subplot(4,1,2)
        spectrogram(data{2}(:,1), hamming(N), N/2, 1024, fs,'yaxis'); % Plots spectrogram
        colormap(jet(fs));  % Choose plot color scheme
        caxis([-80 -40]);   % Choose colomap range (values outside range are clipped to limits)
        title(['Velocidade constante (', speed{2}, ' km/h)']); % Subplot title
        
        % Constant speed 3
        subplot(4,1,3)
        spectrogram(data{3}(:,1), hamming(N), N/2, 1024, fs,'yaxis'); % Plots spectrogram
        colormap(jet(fs));  % Choose plot color scheme
        caxis([-80 -40]);   % Choose colomap range (values outside range are clipped to limits)
        title(['Velocidade constante (', speed{3}, ' km/h)']); % Subplot title
                
        % Accelerating
        subplot(4,1,4)
        spectrogram(data{4}(:,1), hamming(N), N/2, 1024, fs,'yaxis')
        colormap(jet(fs));
        caxis([-80 -40]);
        title('Aceleração');
        
        titleHandle = suptitle(['Espectrograma. Velocidade Constante vs Aceleração. Carro ', carID{1},', sensor 1.']); % Image title
        set(titleHandle,'FontSize',16,'FontWeight','bold')
        F = getframe(gcf);


        case 4
        %% Plot spectrogram: constant speed (all speeds) vs accelerating (cropped)       
        figure;
        set(gcf, 'Position', get(0, 'Screensize'))
        
        [~,~,~,P] = spectrogram(data{1}(:,1), hamming(N), N/2, 1024, fs,'yaxis'); % Signal PSD
        specMin = floor(10*log10(min(P(:))));   % Min and max PSD values to 
        specMax = ceil(10*log10(max(P(:))));    % limit colormap range
        
        % Constant speed 1
        subplot(4,1,1)
        spectrogram(data{1}(:,1), hamming(N), N/2, 1024, fs,'yaxis'); % Plots spectrogram
        colormap(jet(fs));  % Choose plot color scheme
        caxis([-80 -40]);   % Choose colomap range (values outside range are clipped to limits)
        title(['Velocidade constante (', speed{1}, ' km/h)']); % Subplot title
        
        % Constant speed 2
        subplot(4,1,2)
        spectrogram(data{2}(:,1), hamming(N), N/2, 1024, fs,'yaxis'); % Plots spectrogram
        colormap(jet(fs));  % Choose plot color scheme
        caxis([-80 -40]);   % Choose colomap range (values outside range are clipped to limits)
        title(['Velocidade constante (', speed{2}, ' km/h)']); % Subplot title
        
        % Constant speed 3
        subplot(4,1,3)
        spectrogram(data{3}(:,1), hamming(N), N/2, 1024, fs,'yaxis'); % Plots spectrogram
        colormap(jet(fs));  % Choose plot color scheme
        caxis([-80 -40]);   % Choose colomap range (values outside range are clipped to limits)
        title(['Velocidade constante (', speed{3}, ' km/h)']); % Subplot title
                
        % Accelerating
        subplot(4,1,4)
        spectrogram(data{4}(:,1), hamming(N), N/2, 1024, fs,'yaxis')
        colormap(jet(fs));
        caxis([-80 -40]);
        title('Aceleração');
        
        titleHandle = suptitle(['Espectrograma. Velocidade Constante vs Aceleração. Carro ', carID{1},', sensor 1.']); % Image title
        set(titleHandle,'FontSize',16,'FontWeight','bold')
        F = getframe(gcf);
 

end


end

