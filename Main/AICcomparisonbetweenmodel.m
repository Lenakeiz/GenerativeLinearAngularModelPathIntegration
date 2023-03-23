%%
function AICcomparisonbetweenmodel(All_AIC, idx1, idx2, config)

    % plot figures
    f = figure('visible','on','Position', [100 100 200 200]);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)    
    
    AICM1 = All_AIC(:,idx1);
    AICM2 = All_AIC(:,idx2);
    scatter(AICM1, AICM2, 30, 'red', 'filled')
    hold on 
    rline = refline(1,0);
    rline.Color = 'k';
    rline.LineStyle = "--";
    rline.LineWidth = 1;  

    xlabel("AIC M"+num2str(idx1));
    ylabel("AIC M"+num2str(idx2));
    [p, h, stats] = signrank(AICM1, AICM2, "tail","right");
    %print p value and statistics for paper use
    p
    stats
    
    if p<0.001
        title("p<0.001")
    else
        title("p="+num2str(p, '%0.2f'))
    end
    
    
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.0 .0 .0], ...
        'YColor'      , [.0 .0 .0], ...
        'XLim'        , [0,80],...
        'XTick'       , [0,80],...
        'YLim'        , [0,80],...
        'YTick'       , [0,80],...   
        'LineWidth'   , 1        );
    
    exportgraphics(f,config.ResultFolder+"/AICCompare_"+num2str(idx1)+"_"+num2str(idx2)+".png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/AICCompare_"+num2str(idx1)+"_"+num2str(idx2)+".pdf",'Resolution',300,'ContentType','vector');

end