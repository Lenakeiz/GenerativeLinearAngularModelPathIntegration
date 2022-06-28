function PlotActualVsRealreturnDistance(AllX, GroupName)

CreateCustomFigure;
ColorPattern;

savefolder = pwd + "/Output/";
resultfolder = savefolder+"ExtraFigures/ActualVsRealReturnDistance/";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

% Extracting no change condition only
AllX_nochange = AllX{1,1};
AllX_nodistalcue = AllX{1,2};
AllX_noopticflow = AllX{1,3};

slope = [];
intercept = [];

hold on;
for i = 1:length(AllX_nochange)

    % One line per participant
    actualDistances = [];
    idealDistances  = [];

    for j = 1:length(AllX_nochange{1,i})
    currActualDistances = sqrt(...
                            sum(...
                            [(AllX_nochange{1,i}{1,j}(4,1) - AllX_nochange{1,i}{1,j}(3,1))^2, (AllX_nochange{1,i}{1,j}(4,2) - AllX_nochange{1,i}{1,j}(3,2))^2] ...
                            ) ...
                          );
    currIdealDistances = sqrt(...
                            sum(...
                            [(AllX_nochange{1,i}{1,j}(3,1))^2, (AllX_nochange{1,i}{1,j}(3,2))^2] ...
                            ) ...
                          );
    actualDistances = [actualDistances;currActualDistances];
    idealDistances  = [idealDistances;currIdealDistances];
    end
    
    for j = 1:length(AllX_nodistalcue{1,i})
    currActualDistances = sqrt(...
                            sum(...
                            [(AllX_nodistalcue{1,i}{1,j}(4,1) - AllX_nodistalcue{1,i}{1,j}(3,1))^2, (AllX_nodistalcue{1,i}{1,j}(4,2) - AllX_nodistalcue{1,i}{1,j}(3,2))^2] ...
                            ) ...
                          );
    currIdealDistances = sqrt(...
                            sum(...
                            [(AllX_nodistalcue{1,i}{1,j}(3,1))^2, (AllX_nodistalcue{1,i}{1,j}(3,2))^2] ...
                            ) ...
                          );
    actualDistances = [actualDistances;currActualDistances];
    idealDistances  = [idealDistances;currIdealDistances];
    end

    for j = 1:length(AllX_noopticflow{1,i})
    currActualDistances = sqrt(...
                            sum(...
                            [(AllX_noopticflow{1,i}{1,j}(4,1) - AllX_noopticflow{1,i}{1,j}(3,1))^2, (AllX_noopticflow{1,i}{1,j}(4,2) - AllX_noopticflow{1,i}{1,j}(3,2))^2] ...
                            ) ...
                          );
    currIdealDistances = sqrt(...
                            sum(...
                            [(AllX_noopticflow{1,i}{1,j}(3,1))^2, (AllX_noopticflow{1,i}{1,j}(3,2))^2] ...
                            ) ...
                          );
    actualDistances = [actualDistances;currActualDistances];
    idealDistances  = [idealDistances;currIdealDistances];
    end

    randColor = [rand(1,1,1), rand(1,1,1), rand(1,1,1)];
    sp = scatter(idealDistances,actualDistances, SizeData=80, Marker="o",...
                 Color=randColor, MarkerFaceColor=randColor, MarkerEdgeColor=randColor.*0.8,...
                 MarkerFaceAlpha = 0.3, MarkerEdgeAlpha = 0.3);
   
    simpleRegression = polyfit(idealDistances,actualDistances,1); 
    f = polyval(simpleRegression,idealDistances);
        
    slope     = [slope;simpleRegression(1)];
    intercept = [intercept;simpleRegression(2)];

    plot(idealDistances,f,'-', Color=[randColor, 0.4], LineWidth=3);    

end

%Calculate line from
x = 0:0.2:5;
y = x.*mean(slope) + mean(intercept);

plot(x,y,"-",LineWidth=5,Color=[0 , 0, 0, 0.2]);
plot(0:0.2:5,0:0.2:5,"--",LineWidth=2,Color=[0, 0, 0, 0.2]);

hold off;

xlim([0 5]);
ylim([0 5]);

axis square;

ax = gca;

ax.FontSize = 20;

xlabel("Ideal Return Length", "FontSize",20);
ylabel("Actual Return Length", "FontSize",20);

titlelbl = [GroupName];
title(titlelbl,FontSize=25);

exportgraphics(fh,config.ResultFolder+GroupName+".png",'Resolution',300);

end