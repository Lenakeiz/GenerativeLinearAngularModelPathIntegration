function PlotShadedSlidingWindowOfYoung(AllYoungParams, config)
%plot the parameter value against the scale

resultfolder = config.ResultFolder;
Model_Name = config.ModelName;

%create storing folder if not exist
savefoldername = resultfolder+"YoungSliding/"+Model_Name+'/';
if ~exist(savefoldername, 'dir')
   mkdir(savefoldername);
end

numConds = size(AllYoungParams,2);
numSlides = size(AllYoungParams{1},2);
numParams = size(AllYoungParams{1}{1},2);

ResultsMean = zeros(numConds, numSlides, numParams);
ResultsSem = zeros(numConds, numSlides, numParams);
for cond=1:numConds
    ParamsCond = AllYoungParams{cond};
    for slide=1:numSlides
        params = ParamsCond{slide};
        mean_params = mean(params,1,'omitnan');
        sem_params = std(params,1,'omitnan')./sqrt(size(params,1));
        ResultsMean(cond, slide, :) = mean_params;
        ResultsSem(cond, slide, :) = sem_params;
    end
end

%plot for each parameter
title_table = ["bG_1", "bG_2", "bG_3", "g_2", "g_3", "k_3", "sigma", "nu"];
ylimall = [[0.1,1];[0.1,1];[0.1,1];[0.5,2];[0.1,1];[0,pi];[0.1,2];[0.1,100]];
for paramidx=1:numParams
    f = figure('visible','off','Position', [100 100 1000 300]);
    ylabel('Estimated value of parameters');

    subplot(1,3,1);
    y = ResultsMean(1,:,paramidx);
    x = 1:numel(y);
    plot(x, y, '-^k', 'LineWidth', 3);
    hold on
    sem = ResultsSem(1,:,paramidx);
    curve1 = y + sem;
    curve2 = y-sem;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, 'r', 'FaceAlpha',0.5);
    xticks(x);
    ylim(ylimall(paramidx,:));
    ylabel('parameter value');
    xlabel('scale of the outbound path');
    title('No change', 'FontSize', 15);

    subplot(1,3,2);
    y = ResultsMean(2,:,paramidx);
    x = 1:numel(y);
    plot(x, y, '-^k', 'LineWidth', 3);
    hold on
    sem = ResultsSem(2,:,paramidx);
    curve1 = y + sem;
    curve2 = y-sem;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, 'r', 'FaceAlpha',0.5);
    ylim(ylimall(paramidx,:));
    ylabel('parameter value');
    xlabel('scale of the outbound path');
    title('No distal cue', 'FontSize', 15)

    subplot(1,3,3);
    y = ResultsMean(3,:,paramidx);
    x = 1:numel(y);
    plot(x, y, '-^k', 'LineWidth', 3);
    hold on
    sem = ResultsSem(3,:,paramidx);
    curve1 = y + sem;
    curve2 = y-sem;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, 'r', 'FaceAlpha',0.5);  
    ylim(ylimall(paramidx,:));
    ylabel('parameter value');
    xlabel('scale of the outbound path');
    title('No optical flow', 'FontSize', 15)
    
    sgtitle("Parameter name: "+ title_table(paramidx));
    exportgraphics(f,savefoldername+title_table(paramidx)+".png",'Resolution',300);
end
end