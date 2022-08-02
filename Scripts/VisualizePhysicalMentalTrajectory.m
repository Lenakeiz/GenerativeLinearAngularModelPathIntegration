%% Visualize physical and mental trajectory 
% create by Zilong, 01/08/2022

%% Preparing the data
VAM_PrepareBaseConfig;

%% Preprocessing the data
VAM_PreprocessData;

%% Preparing the data and Slecting the Model

config.ModelName        =   "beta_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

% Run the model
VAM;

%% Model run completed, preparing the data for plotting figures
config.ResultFolder     =   pwd + "/Output/ModelFigures/"+config.ModelName+"/VisualizingMentalPhysicalTrajectory";
%create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%% Getting Information from results:
AllYoungResults         =   YoungControls.Results;
AllHealthyOldResults    =   HealthyControls.Results;
AllMCIPosResults        =   MCIPos.Results;
AllMCINegResults        =   MCINeg.Results;
AllMCIUnkResults        =   MCIUnk.Results;

%%
%VisualizeMenPhyTraj(AllHealthyOldResults, 8, 1, 9)
VisualizeMenPhyTraj(AllMCINegResults, 3, 1, 1)
%% Visualize mental physical trajectory
function VisualizeMenPhyTraj(GroupResults, ID, Cond, TrialIdx)

    %extract the parameters
    parameters = GroupResults.estimatedParams{Cond}(ID,:);
    cell_params = num2cell(parameters);
    [beta, g2, g3, sigma, nu] = deal(cell_params{:});

    if isnan(beta)
        error("NAN!")
    end 

    %% extract X
    X = GroupResults.X{Cond}{ID}{TrialIdx}; %X here is sufficient to plot the physical trajectory 
    
    %% extract DX, i.e., leg length
    DX              =       GroupResults.DX{Cond}{ID}{TrialIdx}; 
    l1              =       DX(1);
    l2              =       DX(2);

    %% extract theta
    Theta           =       GroupResults.THETADX{Cond}{ID}{TrialIdx};
    theta2          =       Theta(2); 
    theta3          =       Theta(3); 

    %% find the correct mean return angle based on all trials 
    sampleSize  =   length(GroupResults.DX{1}{ID})+length(GroupResults.DX{2}{ID})+length(GroupResults.DX{3}{ID});
    Alphas      =   zeros(sampleSize,1);
    index       =   0;
    for condddd = 1:3
        trial_num = length(GroupResults.DX{Cond}{ID});
        for trial_id = 1:trial_num
            index = index+1;
            lll_1 = GroupResults.DX{condddd}{ID}{trial_id}(1);
            lll_2 = GroupResults.DX{condddd}{ID}{trial_id}(2);
            thetaaa_2 = GroupResults.THETADX{condddd}{ID}{trial_id}(2);
            
            %calculate the correct return angle
            phy_p1  = [lll_1,0];
            phy_p2  = [lll_1+lll_2*cos(thetaaa_2),lll_2*sin(thetaaa_2)];
            vec1    = phy_p2-phy_p1; vec2 = [0,0]-phy_p2;
            alpha   = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
            alpha   = deg2rad(alpha);%transfer from degree to radians
            alpha   = mod(alpha, 2*pi);  %wrap to (0,2pi)  
            Alphas(index) = alpha;           
        end
    end
    mean_angle = mean(Alphas);

    %% extract duration
    durationL1      =       GroupResults.L1Dur{Cond}{ID}{TrialIdx};
    durationL2      =       GroupResults.L2Dur{Cond}{ID}{TrialIdx};

    %% calculate the mental trajectory
    men_length1     =       l1*(1-exp(-beta*durationL1))/(beta*durationL1)*exp(-beta*durationL2);
    men_p1          =       [men_length1,0];
    
    theta2_prime    =       g2*theta2;

    men_length2     =       l2*(1-exp(-beta*durationL2))/(beta*durationL2);
    men_p2          =       [men_length1+men_length2*cos(theta2_prime),men_length2*sin(theta2_prime)];

    %calculate length of mental vector 3
    h               =       norm(men_p2);
    %calculate turn angle of mental vector 3
    vec1            =       men_p2-men_p1; 
    vec2            =       [0,0]-men_p2;
    alpha           =       atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    alpha           =       deg2rad(alpha);   %transfer from degree to radians
    
    %mental return angle
    sign_alpha      =       sign(alpha);
    theta3_prime    =       g3*abs(alpha)+mean_angle*(1-g3); %reress to mean correct return angle
    theta3_prime    =       sign_alpha*theta3_prime;
        
    
    x3 = h*cos(theta2_prime)*cos(theta3_prime) - h*sin(theta2_prime)*sin(theta3_prime);
    y3 = h*cos(theta2_prime)*sin(theta3_prime) + h*sin(theta2_prime)*cos(theta3_prime);
    men_p3          =       [men_p2(1)+x3, men_p2(2)+y3];


    %% set figure info
    f = figure('visible','on','Position', [100 100 1000 500]);
    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12) 
    
    %plot physical trajectory
    physicalxy  =       X;
    p_x         =       physicalxy(:,1);
    p_y         =       physicalxy(:,2);
    plot(p_x', p_y', '.-', 'Markersize', 10,  LineWidth=1, Color='b');
    hold on

    %plot mental trajectory
    mentalxy    =       [[0,0];men_p1;men_p2;men_p3];
    m_x         =       mentalxy(:,1); 
    m_y         =       mentalxy(:,2);
    plot(m_x', m_y', '.-', 'Markersize', 10,  LineWidth=1, Color='r');

    set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'XLim'        , [-4, 4],...
    'YLim'        , [-4, 4],... 
    'XTick'       , [-4,0,4],...
    'YTick'       , [-4,0,4],...
    'LineWidth'   , .5        );
    xlabel('X (meters)');
    ylabel('Y (meters)');
end