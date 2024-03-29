%% Color Pattern
% Andrea Castegnaro, UCL, 2022 uceeaca@ucl.ac.uk
% This scripts create our color palette and save it in the Output folder
% for visualization. We tried to use specific colors for specific groups
% when possible as outlined below. Storing the information in config
% struct.

% Our color palette
config.color_scheme_npg = [0.9020    0.2941    0.2078; ...
    0.3020    0.7333    0.8353; ...
    0         0.6275    0.5294; ...
    0.2353    0.3294    0.5333; ...
    0.9529    0.6078    0.4980; ...
    0.5176    0.5686    0.7059; ...
    0.5686    0.8196    0.7608; ...
    0.8627         0         0; ...
    0.4941    0.3804    0.2824; ...
    0.6902    0.6118    0.5216 ];

% In the following order: Young, old, mci unk(grouped), mci neg, mci
% pos
config.color_scheme_group = config.color_scheme_npg([1 5 2 4 6],:);

folder = pwd + "/Output/ColorPattern";
%create storing folder for trajectory if not exist
if ~exist(folder, 'dir')
    mkdir(folder);
end

f = figure('Visible','off','Position',get(0,'ScreenSize'));

hold on;
for i = 1:length(config.color_scheme_npg)
    rectangle('Position',[(i-1)*2 0 2 2],'FaceColor',config.color_scheme_npg(i,:));
end

text(1,1,'Young', 'Rotation',90, 'FontSize',25);
text(3,1,'MCI Grouped or unknown', 'Rotation',90, 'FontSize',25);
text(7,1,'MCI negative', 'Rotation',90, 'FontSize',25);
text(11,1,'MCI positive', 'Rotation',90, 'FontSize',25);
text(9,1,'Healthy Older', 'Rotation',90, 'FontSize',25);

text(15,0.1,'Gamma', 'Rotation',90, 'FontSize',25);
text(5,0.1,'g', 'Rotation',90, 'FontSize',25);
text(13,0.1,'k', 'Rotation',90, 'FontSize',25);
text(17,0.1,'sigma', 'Rotation',90, 'FontSize',25);
text(19,0.1,'niu', 'Rotation',90, 'FontSize',25);
hold off;

title("Color Scheme", "FontSize",30)
axis([0 length(config.color_scheme_npg)*2 0 2])
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

exportgraphics(f,folder+"/ColorPattern.png",'Resolution',300);

clear folder f i