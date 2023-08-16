%% Script to visualize the distribution of l1 l2 and theata2 in different groups
% Zilong Ji, UCL, 2022 zilong.ji@ucl.ac.uk

% Preparing the data
GLAMPI_PrepareBaseConfig

% Preprocessing the data
GLAMPI_PreprocessData

%%

% Model fitting
config.ModelName        =   "beta_k_g2_g3_m3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "m3", "sigma", "nu"];

config.NumParams        =   length(config.ParamName);
GLAMPI;

%%
% Generating color scheme for our paper
ColorPattern; 

%%
%get l1 l2 and theta2
[YoungL1, YoungL2] = extractLengthData(YoungControls);
[HOL1, HOL2] = extractLengthData(HealthyControls);
[MCIUnkL1, MCIUnkL2] = extractLengthData(MCIUnk);
[MCIPosL1, MCIPosL2] = extractLengthData(MCIPos);
[MCINegL1, MCINegL2] = extractLengthData(MCINeg);

YoungA2 = extractAngleData(YoungControls);
HOA2 = extractAngleData(HealthyControls);
MCIUnkA2 = extractAngleData(MCIUnk);
MCIPosA2 = extractAngleData(MCIPos);
MCINegA2 = extractAngleData(MCINeg);


config.ResultFolder = pwd + "/Output/FigS10";
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%%
histgramplot(YoungL1, HOL1, MCIUnkL1, MCIPosL1, MCINegL1, "Leg 1", config)
histgramplot(YoungL2, HOL2, MCIUnkL2, MCIPosL2, MCINegL2, "Leg 2", config)
histgramplot(YoungA2, HOA2, MCIUnkA2, MCIPosA2, MCINegA2, "Angle", config)

%%
function histgramplot(Young, HO, MCIUnk, MCIPos, MCINeg, type, config)

    f = figure('visible','off','Position', [100 100 500 700]);

    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)     

    colorY = config.color_scheme_npg(1,:);
    colorH = config.color_scheme_npg(2,:);
    colorU = config.color_scheme_npg(3,:);
    colorP = config.color_scheme_npg(4,:);
    colorN = config.color_scheme_npg(5,:);

    subplot(5,1,1)
    histogram(Young, 50, 'FaceColor', colorY);
    title('Young');
    box off; 
    if strcmp(type, "Leg 1") || strcmp(type, "Leg 2")
        xlim([1,5]);
    else
        xlim([0,3]);
    end

    subplot(5,1,2)
    histogram(HO, 50, 'FaceColor', colorH);
    title('Healthy Elderly');
    box off; 
    if strcmp(type, "Leg 1") || strcmp(type, "Leg 2")
        xlim([1,5]);
    else
        xlim([0,3]);
    end

    subplot(5,1,3)
    histogram(MCIUnk, 50, 'FaceColor', colorU);
    title('MCIUnk');
    box off; 
    if strcmp(type, "Leg 1") || strcmp(type, "Leg 2")
        xlim([1,5]);
    else
        xlim([0,3]);
    end

    subplot(5,1,4)
    histogram(MCIPos, 50, 'FaceColor', colorP);
    title('MCIPos');
    box off; 
    if strcmp(type, "Leg 1") || strcmp(type, "Leg 2")
        xlim([1,5]);
    else
        xlim([0,3]);
    end

    subplot(5,1,5)
    histogram(MCINeg, 50, 'FaceColor', colorN);
    title('MCINeg');
    box off; 
    if strcmp(type, "Leg 1") || strcmp(type, "Leg 2")
        xlim([1,5]);
    else
        xlim([0,3]);
    end
    
    if strcmp(type, "Leg 1") || strcmp(type, "Leg 2")
        xlabel(type+" (meters)");
    else
        xlabel(type+" (radians)");
    end
    
    exportgraphics(f,config.ResultFolder+"/"+type+".png",'Resolution',300);
end

%%
function [L1_all, L2_all] = extractLengthData(groupData)

    L1_all = []; L2_all = [];
    for cond=1:3
        L1_cond = []; L2_cond = [];
        for i_participant = 1:length(groupData.Results.DX{cond})
            if(~isempty(groupData.Results.DX{cond}{i_participant}))
                paths = cell2mat(groupData.Results.DX{cond}{i_participant});
                %extract path 1 and path 2
                l1 = paths(1,:);
                l2 = paths(2,:);
                L1_cond = [L1_cond, l1];
                L2_cond = [L2_cond, l2];
            end
        end
        L1_all = [L1_all, L1_cond];
        L2_all = [L2_all, L2_cond];
    end
end

function A2_all = extractAngleData(groupData)
    
    A2_all = [];
    for cond=1:3
        A2_cond = [];
        for i_participant = 1:length(groupData.Results.THETADX{cond})
            if(~isempty(groupData.Results.THETADX{cond}{i_participant}))
                paths = cell2mat(groupData.Results.THETADX{cond}{i_participant});
                %extract path 1 and path 2
                a2 = paths(2,:);
                A2_cond = [A2_cond, a2];
            end
        end
        A2_all = [A2_all, A2_cond];
    end
end
