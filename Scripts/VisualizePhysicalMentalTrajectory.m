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

%% Visualize mental physical trajectory
function VisualizeMenPhyTraj(GroupResults, ID, Cond, TrialIdx)
    %extract the parameters
    parameters = GroupResults.estimatedParams{Cond}(ID,:);
    cell_params = num2cell(parameters);
    [beta, g2, g3, sigma, nu] = deal(cell_params{:});

    %extract X
    X = GroupResults.X{Cond}{ID}{TrialIdx}; %X here is sufficient to plot the physical trajectory 
    
    %extract DX, i.e., leg length
    DX = GroupResults.DX{Cond}{ID}{TrialIdx}; 
    l1              =       DX(1);
    l2              =       DX(2);

    %extract theta
    Theta           =       GroupResults.THETADX{Cond}{ID}{TrialIdx};
    theta2          =       THETAX(2); 
    theta3          =       THETAX(3); 

    %extract duration
    L1Dur           =       GroupResults.L1Dur{Cond}{ID}{TrialIdx};
    L2Dur           =       GroupResults.L2Dur{Cond}{ID}{TrialIdx};

    %calculate the mental trajectory
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

    men_p3          =       [men_p2(1)+h*cos(theta3_prime), men_p2(2)+h*sin(theta3_prime)];


    %% set figure info
    f = figure('visible','off','Position', [100 100 1000 500]);
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
end