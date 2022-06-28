%% Two-way anova on real data: distance, angle or mean walking speed
% five groups:         Young / HealthyOld / MCIPos / MCINeg / MCIUnk
% three conditions:    no change / no distal cue /  no optical flow
function ThreewayAnova_CocoData(PosData, NegData, config)

    %load configurations necessary for the script
    resultfolder = config.ResultFolder;

    %create storing folder for trajectory if not exist
    savefoldername = resultfolder+"/ThreewayAnova/";
    if ~exist(savefoldername, 'dir')
       mkdir(savefoldername);
    end

    PosGender = config.Gender_Pos;
    NegGender = config.Gender_Neg;

    %calculate the mean error for each participant
    %swap last two condition
    Condition = [1,3,2];
    
    PosMean = cell(1,3);
    for i=1:3
        cond = Condition(i);
        PosMean_cond_i = PosData{cond};
        Mean = [];
        for id=1:length(PosMean_cond_i)
            err = PosMean_cond_i{id};
            err = cell2mat(err);
            meanerr = mean(err, 'omitnan');
            Mean = [Mean, meanerr];
        end
        PosMean{i} = Mean;
    end

    NegMean = cell(1,3);
    for i=1:3
        cond = Condition(i);
        NegMean_cond_i = NegData{cond};
        Mean = [];
        for id=1:length(NegMean_cond_i)
            err = NegMean_cond_i{id};
            err = cell2mat(err);
            meanerr = mean(err, 'omitnan');
            Mean = [Mean, meanerr];
        end
        NegMean{i} = Mean;
    end

    %processing the data into a long numeric vector 
    [PosY, PosGroupNames, PosConditionNames, PosGenderNames]=ReGroupDataForAnova(PosMean,'Pos', PosGender);
    [NegY, NegGroupNames, NegConditionNames, NegGenderNames]=ReGroupDataForAnova(NegMean,'Neg', NegGender);

    AllY = [PosY,NegY];
    AllGroupNames = [PosGroupNames,NegGroupNames];
    AllConditionNames = [PosConditionNames,NegConditionNames];
    AllGenderNames = [PosGenderNames, NegGenderNames];

    %Do three-way anova with unbalanced design
    [p,anova_tab, stats]= anovan(AllY,{AllGroupNames,AllConditionNames, AllGenderNames},'model','full','varnames',{'Groups','Conditions','Genders'},'display','on');
    
    %Do multiple comparisons on main effect 1
    figure()
    multcompare(stats,'Dimension',[1],'CType','bonferroni');
    saveas(gcf,savefoldername+"ME1.png");
    close(gcf);

    %Do multiple comparisons on main effect 2
    figure()
    multcompare(stats,'Dimension',[2],'CType','bonferroni');
    saveas(gcf,savefoldername+"ME2.png");
    close(gcf);

    %Do multiple comparisons on main effect 3
    figure()
    multcompare(stats,'Dimension',[3],'CType','bonferroni');
    saveas(gcf,savefoldername+"ME3.png");
    close(gcf);

    %Do multiple comparisons on main effect 1&2
    figure()
    multcompare(stats,'Dimension',[1 2],'CType','bonferroni');
    saveas(gcf,savefoldername+"ME1ME2.png");
    close(gcf);    

    %Do multiple comparisons on main effect 1&3
    figure()
    multcompare(stats,'Dimension',[1 3],'CType','bonferroni');
    saveas(gcf,savefoldername+"ME1ME3.png");
    close(gcf); 

    %Do multiple comparisons on main effect 2&3
    figure()
    multcompare(stats,'Dimension',[2 3],'CType','bonferroni');
    saveas(gcf,savefoldername+"ME2ME3.png");
    close(gcf); 

    %Do multiple comparisons on main effect 1&2&3
    figure()
    multcompare(stats,'Dimension',[1 2 3],'CType','bonferroni');
    saveas(gcf,savefoldername+"ME1ME2ME3.png");
    close(gcf);    
end

function [Y, GroupNames, ConditionNames, GenderNames] = ReGroupDataForAnova(Data, groupname, gender)
    Y = [Data{1}, Data{2}, Data{3}];

    numSubjs = size(Data{1},2);
    
    GroupNames = [string(repmat({groupname},1,numSubjs)),...
                  string(repmat({groupname},1,numSubjs)),...
                  string(repmat({groupname},1,numSubjs))];
    
    ConditionNames = [string(repmat({'NoChange'},1,numSubjs)),...
                      string(repmat({'NoOpticFlow'},1,numSubjs)),...
                      string(repmat({'NoDistalCue'},1,numSubjs))];  

    GenderNames = [gender', gender', gender']; 
end