function fig = plot_Cmat(t, tau, Cmat, salvar, filename, titulo)
    figurebackcolor = 'white';
%     pos = [0.01 0.52 0.49 0.40];
    fp1 = figure('numbertitle','off','name','GCC e TDD',...
             'Units','normal');
    colordef(fp1,figurebackcolor);
    
    % Plot correlation matrix as image
    imagesc(t, 1000*tau, Cmat); 
    colormap(flipud(bone)); colorbar; 
    set(gca,'YDir','normal'); hold on
    set(gcf,'Visible', 'off');
    
    % Edit figure
    hAx = gca; 
    ylabel(gca, 'Delay \tau (ms)') % left y-axis
    xlabel('Time (s)');
    title(titulo, 'Interpreter', 'none');
    grid on;  
%     p = [hLines(1); hLines(3)];
%     legenda = {'Estimated TDoA', 'Theoretical TDoA'};
%     legend('show'); legend(p, legenda,'Location','best');
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