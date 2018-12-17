function fig = plot_tdd(t, tau, tdd, Cmat, salvar, filename, titulo)
    figurebackcolor = 'white';
    pos = [0.01 0.52 0.49 0.40];
    fp1 = figure('numbertitle','off','name','GCC e TDD',...
             'Units','normal','Position', pos);
    colordef(fp1,figurebackcolor);
    
    % Plot correlation matrix as image
    imagesc(t,1000*tau,Cmat); 
    colormap(flipud(bone)); %colorbar; 
    set(gca,'YDir','normal'); hold on
    set(gcf,'Visible', 'off');
    
    % Plot tdd curve
    p = plot(t, 1000*tdd, 'Color','blue', 'Linewidth', 1);
    ylabel('Atraso \tau (ms)');
    xlabel('Tempo (s)');
    
    % Edit figure
    title(titulo);
    %yticks([-0.6:0.2:0.6])
    grid on;
    legend('show'); legend(p,   'TDD Estimado', ...
                                'Location', 'best');
    set(gca,'YTick',[-0.8:0.2:0.8])
    fig = gcf;
    fig.InvertHardcopy = 'off';
    
    if salvar == true
%        print(gcf, filename,'-dpng','-r300'); % Salvando imagem em arquivo png
        print(gcf, filename,'-depsc','-r300'); % Salvando imagem em arquivo eps
    end
end