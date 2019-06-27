function fig = plot_theo_compare(t, tau, tdd, salvar, filename, titulo)
    figurebackcolor = 'white';
%     pos = [0.01 0.52 0.49 0.40];
    fp1 = figure('numbertitle','off','name','GCC e TDD',...
             'Units','normal');
    colordef(fp1,figurebackcolor);
    
%     % Plot correlation matrix as image
%     imagesc(t, 1000*tau, Cmat); 
%     colormap(flipud(bone)); %colorbar; 
%     set(gca,'YDir','normal'); hold on
    set(gcf,'Visible', 'off');
    
    % Plot tdd curve
    hLines = plot([t', t', t', t'], 1000*tdd);
    hAx = gca; 
    ylabel(gca, 'Delay \tau (ms)') % left y-axis
    xlabel('Time (s)');
    hLines(1).Color = 'blue';
    hLines(2).Color = 'blue';
    hLines(3).Color = 'red';
    hLines(4).Color = 'red';
    hLines(1).LineStyle = ':';
    hLines(2).LineStyle = ':';
    hLines(3).LineStyle = ':';
    hLines(4).LineStyle = ':';
    hLines(1).LineWidth = 2;
    hLines(2).LineWidth = 2;
    hLines(3).LineWidth = 2;
    hLines(4).LineWidth = 2;
    
    % Edit figure
    title(titulo, 'Interpreter', 'none');
    grid on;  
    p = [hLines(1); hLines(3)];
    legenda = {'Mean distance', 'Measured distance'};
    legend('show'); legend(p, legenda,'Location','best');
    set(hAx,'YTick',[-0.8:0.1:0.8])
    maxTau = 1000*max(abs(tau));
    set(hAx,'ylim',[-maxTau maxTau])
    set(hAx,'xlim',[0 t(end)])
    set(hAx,'ycolor','k')
    fig = gcf;
    fig.InvertHardcopy = 'off';
    
    if salvar == true
%        print(gcf, filename,'-dpng','-r300'); % Salvando imagem em arquivo png
        print(gcf, filename,'-depsc','-r300'); % Salvando imagem em arquivo eps
    end
end