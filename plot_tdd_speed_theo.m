function fig = plot_tdd_speed_theo(t, tau, tdd, Cmat, salvar, filename, titulo, speed)
    figurebackcolor = 'white';
    pos = [0.01 0.52 0.49 0.40];
    fp1 = figure('numbertitle','off','name','GCC e TDD',...
             'Units','normal','Position', pos);
    colordef(fp1,figurebackcolor);
%     left_color = RGB('black');
%     right_color = RGB('purple');
%     set(fp1,'defaultAxesColorOrder',[left_color; right_color]);
    
    % Plot correlation matrix as image
    imagesc(t, 1000*tau, Cmat); 
    colormap(flipud(bone)); %colorbar; 
    set(gca,'YDir','normal'); hold on
    set(gcf,'Visible', 'off');
            
    % Plot tdd curve
    [hAx, hLeft, hRigth] = plotyy([t', t', t', t'], 1000*tdd, t, speed);
    
    ylabel(hAx(1),'Delay \tau (ms)') % left y-axis
    ylabel(hAx(2),'Speed (km/h)') % right y-axis
    xlabel('Tempo (s)');
    hLeft(1).Color = 'blue';
    hLeft(2).Color = 'blue';
    hLeft(3).Color = 'red';%RGB('purple');
    hLeft(4).Color = 'red';%RGB('purple');
    hLeft(1).LineWidth = 2;
    hLeft(2).LineWidth = 2;
    hLeft(3).LineWidth = 2;
    hLeft(4).LineWidth = 2;
    hRigth(1).LineWidth = 2;
    hRigth(1).Color = RGB('purple');
    hLeft(1).LineStyle = '--';
    hLeft(2).LineStyle = '--';
    hLeft(3).LineStyle = '-.';
    hLeft(4).LineStyle = '-.';
        
    % Edit figure
    title(titulo, 'Interpreter', 'none');
    grid on;  
    p = [hLeft(1); hLeft(3); hRigth];
    legenda = {'Estimated TDoA', 'Theoretical TDoA', 'Speed'};
    legend('show'); legend(p, legenda,'Location','best');
    maxTau = 1000*max(abs(tau));
    set(hAx(1),'YTick',[-0.8:0.2:0.8])
    set(hAx(2),'YTick',[0:20:100])
    set(hAx(1),'ylim',[-maxTau maxTau])
    set(hAx(2),'ylim',[0 100])
    set(hAx(1),'xlim',[0 t(end)])
    set(hAx(2),'xlim',[0 t(end)])
    set(hAx(1),'ycolor','black')
    set(hAx(2),'ycolor', RGB('purple'))
    fig = gcf;
    fig.InvertHardcopy = 'off';
    
    if salvar == true
%        print(gcf, filename,'-dpng','-r300'); % Salvando imagem em arquivo png
        print(gcf, filename,'-depsc','-r300'); % Salvando imagem em arquivo eps
    end
end