%% ------------------------- Reset the system ----------------------------%%
close all

clear all
clc

tic

%Time Unit: hour
dt=1;
tmin=0;
tmax=1000;
tspan = [tmin:dt:tmax];

%EDC concentrations
X_ini = 0.01; %Initial non-zero concentration
X_final = 1E6; %Final concentration
X_interval = 1.05; %Geometric increment
%

%% ------------------------- PARAMETERS ----------------------------------%%
load('Default_param.mat');
param = default_param;


%% ------------------------- INITIAL CONDITION----------------------------%%

init_EH = 10;
init_PH = 1;
init_EHCR = 0;
init_XCR = 0;
init_EHPR = 0;
init_XPR = 0;

y0 = [init_EH, init_PH, init_EHCR, init_XCR, init_EHPR, init_XPR];

%% ------------------------- TIME-COURSE SIMULATION (Fig. S1) ------------%%
param = default_param;
y0 = [init_EH, init_PH, init_EHCR, init_XCR, init_EHPR, init_XPR];

[t,y] = ode23s('EDC_ode',tspan, y0, [], param); 
y0 = y(end,:);

tspan1 = [0:0.5:25];
[t1,y1] = ode23s('EDC_ode',tspan1, y0, [], param); 
y0 = y1(end,:);

% For primary hyper- or hypo-functioning, change k1 to 5-fold above or below its default value; for secondary hyper- or hypo-functioning, change k3 to 5-fold above or below its default value.
param.k1 = 5*default_param.k1;
tspan2 = [25:0.5:100];
[t2,y2] = ode23s('EDC_ode',tspan2, y0, [], param); 


%Plotting EH
figure(100)
plot([t1;t2], [y1(:,1); y2(:,1)], 'LineWidth', 3)
hold on
xlim([0,100])
ylim([0,20])
xlabel('Time (h)')
ylabel('EH')
set(gca,'FontSize',18)
pbaspect([1.8 1 1])
box off

%Plotting PH"
figure(101)
plot([t1;t2], [y1(:,2); y2(:,2)], 'LineWidth', 3)
hold on
xlim([0,100])
ylim([0,6])
yticks([0:1:6])
xlabel('Time (h)')
ylabel('PH')
set(gca,'FontSize',18)
pbaspect([1.8 1 1])
box off



%% ------------------------- SINGLE-PARAMETER SIMULATION FOR X AS AN AGONIST (Figs. 2-4, Figs. S3-S4) --------------------------------%%

y0 = [init_EH, init_PH, init_EHCR, init_XCR, init_EHPR, init_XPR];
param = default_param;

%Parameters to vary: select one parameter at a time from the list below and then run this entire section
param_name = "kf6";
% kf6 (kb6) - Fig. 2A-2C, 3A
% kf5 (kb5) - Fig. 2D-2F, 3B
% kf8 (kb8) - Fig. 2G-2I, 3C
% kf7 (kb7) - Fig. 2J-2L, 3D
% wp        - Fig. 4A-4C, 3E
% wc        - Fig. 4D-4F, 3F
% k1, k2, k3, k4, k30, Kd3, n3, CRtot, or PRtot - Monotonic dose-response results
% none      - Fig. S3

%Range of fold change of the parameter
if param_name == "none" %For Fig. S3
    fold_change_vec = [1];
else
    if param_name == "n3" %For Fig. S4
        fold_change_vec = [1000/7, 100/7, 10/7, 1, 4/7, 1/7];
    else
        fold_change_vec = [8,4,2,1,0.5,0.25,0.125];
    end

    %Actual values of the parameter for the above fold changes
    param_value_vec = eval(strcat('default_param.', param_name, '*fold_change_vec'));
    eval(strcat(param_name, '_vec = [', num2str(param_value_vec), '];'));

end

%Number of fold changes
P_length = length(fold_change_vec);


%To set basal CR saturation at 1% while keeping the basal PH and EH at same levels as when the saturation is the default 10%
% param.kb7 = 990;
% param.CRtot = 100;

%To set basal CR saturation at 50% while keeping the basal PH and EH at same levels as when the saturation is the default 10%
% param.kb7 = 10;
% param.CRtot = 2;


%Defining color shade for plotting
bluecolor = '#0072BD';%'#004F83';%'#0072BD' ;
nsteps = P_length + 2;
% convert to a usable RGB tuple
color = sscanf(bluecolor(2:end),'%2x%2x%2x',[1 3])/255;
% calculate output colors
scale = linspace(1,0,nsteps);
blueshades = scale .* color(:);
color_inv = (1-color);
bluetints =  fliplr(scale) .* color_inv(:) + color(:);

redcolor = '#FF0000';%'#D95319' ;
nsteps = P_length + 2;
% convert to a usable RGB tuple
color = sscanf(redcolor(2:end),'%2x%2x%2x',[1 3])/255;
% calculate output colors
scale = linspace(1,0,nsteps);
redshades = scale .* color(:);
color_inv = (1-color);
redtints =  fliplr(scale) .* color_inv(:) + color(:);

greencolor = '#77AA30';%'#77AA30' ;
nsteps = P_length + 2;
% convert to a usable RGB tuple
color = sscanf(greencolor(2:end),'%2x%2x%2x',[1 3])/255;
% calculate output colors
scale = linspace(1,0,nsteps);
greenshades = scale .* color(:);
color_inv = (1-color);
greentints =  fliplr(scale) .* color_inv(:) + color(:);

%Looping through different parameter values
for j = 1:1:P_length

    j

    if param_name ~= "none"
        eval(strcat('param.', param_name, '=', num2str(param_value_vec(j)), ';'));
    end

    %For Fig. S4
    if param_name == "n3" 
        param.Kd3 = (0.09/(1-0.09))^(1/param.n3); %for the case of varying n3 while keeping k3-mediated rate at 9% of k3.
    end
    
    %EDC concentration
    param.X = 0;
    
    i = 1;
    while (param.X <= X_final)  

        options2 = odeset('RelTol',1e-8,'AbsTol',1e-20);
        [t,y] = ode23tb('EDC_ode',tspan, y0, options2, param); 

        %Steay-state values of model output
        EHss(i)=y(end,1);
        PHss(i)=y(end,2);
        EHCRss(i)=y(end,3);
        XCRss(i)=y(end,4);
        EHPRss(i)=y(end,5);
        XPRss(i)=y(end,6);
        CRss(i)=param.CRtot-EHCRss(i)-XCRss(i);
        PRss(i)=param.PRtot-EHPRss(i)-XPRss(i);
        Endocrine_Effect(i)=EHPRss(i)+param.wp*XPRss(i);
        Endocrine_Effect_BasalValue_vec(j) = Endocrine_Effect(1);
        Normalized_Endocrine_Effect(i) = Endocrine_Effect(i)/Endocrine_Effect(1);
        X_vec(i) = param.X;
        
        
        % %Time course to observe whether steady state is reached
        % figure(100)
        % plot(t,y(:,1))
        % hold on
        % 
        % figure(200)
        % plot(t,y(:,6))
        % hold on

        %Increase X concentration
        if (param.X == 0)
            param.X = X_ini;
        else
            param.X = param.X * X_interval;
        end
        
        i = i + 1

    end %End of while loop

      
    %PH
    figure(1)
    pbaspect([1 1 1])
    semilogx(X_vec,PHss,'LineWidth',3)
    xlabel('Agonist X)','FontSize',16,'FontWeight','bold')
    ylabel('PH','FontSize',16,'FontWeight','bold')
    xlim([0.01 100000])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northeast')
    legend boxon  

    %EH
    figure(2)
    pbaspect([1 1 1])
    semilogx(X_vec,EHss,'LineWidth',3)
    xlabel('Agonist X','FontSize',16,'FontWeight','bold')
    ylabel('EH','FontSize',16,'FontWeight','bold')
    xlim([0.01 100000])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northeast')
    legend boxon  
    
    %EHCR
    figure(3)
    pbaspect([1 1 1])
    semilogx(X_vec,EHCRss,'LineWidth',3)
    xlabel('Agonist X','FontSize',16,'FontWeight','bold')
    ylabel('EHCR','FontSize',16,'FontWeight','bold')
    xlim([0.01 100000])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northeast')
    legend boxon 

    %XCR
    figure(4)
    pbaspect([1 1 1])
    semilogx(X_vec,XCRss,'LineWidth',3)
    xlabel('Agonist X','FontSize',16,'FontWeight','bold')
    ylabel('XCR','FontSize',16,'FontWeight','bold')
    xlim([0.01 100000])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northwest')
    legend boxon

    %EHPR
    figure(5)
    pbaspect([1 1 1])
    semilogx(X_vec,EHPRss,'LineWidth',3)
    xlabel('Agonist X')
    ylabel('EHPR')
    xlim([0.01 100000])
    xticks([1E-2 1E0 1E2 1E4 1E6])
    set(gca,'FontSize',22)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northeast','FontSize',17)
    legend boxon

    %XPR
    figure(6)
    pbaspect([1 1 1])
    semilogx(X_vec,XPRss,'LineWidth',3)
    xlabel('Agonist X')
    ylabel('XPR')
    xlim([0.01 100000])
    xticks([1E-2 1E0 1E2 1E4 1E6])
    set(gca,'FontSize',22)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northwest','FontSize',17)
    legend boxon
    
    %Endocrine Effect
    figure(7)
    pbaspect([2.2 1 1])
    semilogx(X_vec,Endocrine_Effect,'LineWidth',3)
    xlabel('Agonist X')
    ylabel('Endocrine Effect')
    xlim([0.01 100000])
    xticks([1E-2 1E-1 1E0 1E1 1E2 1E3 1E4 1E5 1E6])
    set(gca,'FontSize',15)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northwest','FontSize',15)
    legend boxon 

    %Normalized Endocrine Effect
    figure(8)
    pbaspect([1 1 1])
    semilogx(X_vec,Normalized_Endocrine_Effect,'LineWidth',3)
    xlabel('Agonist X','FontSize',16,'FontWeight','bold')
    ylabel('Normalized Endocrine Effect','FontSize',16,'FontWeight','bold')
    xlim([0.01 100000])
    xticks([1E-2 1E-1 1E0 1E1 1E2 1E3 1E4 1E5 1E6])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northwest')
    legend boxon 

    %EE,  wp*XPR, and EHPR
    figure(9)
    pbaspect([1.5 1 1])
    semilogx(X_vec,Endocrine_Effect,'LineWidth',4,'color',bluetints(:,j))
    hold on
    %A different legend name, legendInfo3, is used to avoid mix-up with the
    %name above because here we are plotting 3 lines at a time.
    legendInfo3{(j-1)*3+1} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    
    semilogx(X_vec,param.wp*XPRss,':','LineWidth',3,'color',redtints(:,j))
    legendInfo3{(j-1)*3+2} = [char(param_name), ' x ', num2str(fold_change_vec(j))];

    semilogx(X_vec,EHPRss,':','LineWidth',3,'color',greentints(:,j))
    legendInfo3{(j-1)*3+3} = [char(param_name), ' x ', num2str(fold_change_vec(j))];

    xlabel('X (agonist)','FontSize',16')
    ylabel('EE, wp*XPR, and EHPR','FontSize',16)
    xlim([0.01 100000])
    xticks([1E-2 1E-1 1E0 1E1 1E2 1E3 1E4 1E5 1E6])
    %set(gca,'FontSize',16,'LineWidth',3)
    set(gca,'FontSize',12)
    box off
    legend(legendInfo3,'location','northeast')
    legend boxon
    
    %CR
    figure(10)
    pbaspect([1 1 1])
    semilogx(X_vec,CRss,'LineWidth',3)
    xlabel('Agonist X','FontSize',16,'FontWeight','bold')
    ylabel('CR','FontSize',16,'FontWeight','bold')
    xlim([0.01 100000])
    xticks([1E-2 1E0 1E2 1E4 1E6])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northwest')
    legend boxon
    
    %PR
    figure(11)
    pbaspect([1 1 1])
    semilogx(X_vec,PRss,'LineWidth',3)
    xlabel('Agonist X','FontSize',16,'FontWeight','bold')
    ylabel('PR','FontSize',16,'FontWeight','bold')
    xlim([0.01 100000])
    xticks([1E-2 1E0 1E2 1E4 1E6])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northwest')
    legend boxon

    
end %End of j for loop


%% ------------------------- SINGLE-PARAMETER SIMULATION FOR X AS AN ANTAGONIST (Figs. 5-6, Figs. S5-S7) -----------------------------%%

y0 = [init_EH, init_PH, init_EHCR, init_XCR, init_EHPR, init_XPR];
param = default_param;

%Set wp and wc to 0 for antagonist simulations
param.wp = 0;
param.wc = 0;

%Parameters to vary: select one parameter at a time from the list below and then run this entire section
param_name = "kf6";
% kf6 (kb6) - Fig. 5A-5C,   6A
% kf5 (kb5) - Fig. S7A-S7C, 6B
% kf8 (kb8) - Fig. 5D-5F,   6C
% kf7 (kb7) - Fig. S7D-S7G, 6D
% k1, k2, k3, k4, k30, Kd3, n3, CRtot, or PRtot - Monotonic dose-response results
% none      - Fig. S5

%Range of fold change of the parameter
if param_name == "none" %For Fig. S5
    fold_change_vec = [1];
else
    if param_name == "n3" %For Fig. S6
        fold_change_vec = [1000/7, 100/7, 10/7, 1, 4/7, 1/7];
    else
        fold_change_vec = [8,4,2,1,0.5,0.25,0.125];
    end

    %Actual values of the parameter for the above fold changes
    param_value_vec = eval(strcat('default_param.', param_name, '*fold_change_vec'));
    eval(strcat(param_name, '_vec = [', num2str(param_value_vec), '];'));

end

%Number of fold changes
P_length = length(fold_change_vec);


%To set basal CR saturation at 1% while keeping the basal PH and EH at same levels as when the saturation is the default 10%
% param.kb7 = 990;
% param.CRtot = 100;

%To set basal CR saturation at 50% while keeping the basal PH and EH at same levels as when the saturation is the default 10%
% param.kb7 = 10;
% param.CRtot = 2;


%Defining color shade for plotting
bluecolor = '#0072BD';%'#004F83';%'#0072BD' ;
nsteps = P_length + 2;
% convert to a usable RGB tuple
color = sscanf(bluecolor(2:end),'%2x%2x%2x',[1 3])/255;
% calculate output colors
scale = linspace(1,0,nsteps);
blueshades = scale .* color(:);
color_inv = (1-color);
bluetints =  fliplr(scale) .* color_inv(:) + color(:);

redcolor = '#FF0000';%'#D95319' ;
nsteps = P_length + 2;
% convert to a usable RGB tuple
color = sscanf(redcolor(2:end),'%2x%2x%2x',[1 3])/255;
% calculate output colors
scale = linspace(1,0,nsteps);
redshades = scale .* color(:);
color_inv = (1-color);
redtints =  fliplr(scale) .* color_inv(:) + color(:);

greencolor = '#77AA30';%'#77AA30' ;
nsteps = P_length + 2;
% convert to a usable RGB tuple
color = sscanf(greencolor(2:end),'%2x%2x%2x',[1 3])/255;
% calculate output colors
scale = linspace(1,0,nsteps);
greenshades = scale .* color(:);
color_inv = (1-color);
greentints =  fliplr(scale) .* color_inv(:) + color(:);

%Looping through different parameter values
for j = 1:1:P_length

    j

    if param_name ~= "none"
        eval(strcat('param.', param_name, '=', num2str(param_value_vec(j)), ';'));
    end

    %For Fig. S6
    if param_name == "n3" 
        param.Kd3 = (0.09/(1-0.09))^(1/param.n3); %for the case of varying n3 while keeping k3-mediated rate at 9% of k3.
    end
    
    %EDC concentration
    param.X = 0;
    
    i = 1;
    while (param.X <= X_final)  

        options2 = odeset('RelTol',1e-8,'AbsTol',1e-20);
        [t,y] = ode23tb('EDC_ode',tspan, y0, options2, param); 

        %Steay-state values of model output
        EHss(i)=y(end,1);
        PHss(i)=y(end,2);
        EHCRss(i)=y(end,3);
        XCRss(i)=y(end,4);
        EHPRss(i)=y(end,5);
        XPRss(i)=y(end,6);
        CRss(i)=param.CRtot-EHCRss(i)-XCRss(i);
        PRss(i)=param.PRtot-EHPRss(i)-XPRss(i);
        Endocrine_Effect(i)=EHPRss(i)+param.wp*XPRss(i);
        Endocrine_Effect_BasalValue_vec(j) = Endocrine_Effect(1);
        Normalized_Endocrine_Effect(i) = Endocrine_Effect(i)/Endocrine_Effect(1);
        X_vec(i) = param.X;
        
        
        % %Time course to observe whether steady state is reached
        % figure(100)
        % plot(t,y(:,1))
        % hold on
        % 
        % figure(200)
        % plot(t,y(:,6))
        % hold on

        %Increase X concentration
        if (param.X == 0)
            param.X = X_ini;
        else
            param.X = param.X * X_interval;
        end
        
        i = i + 1

    end %End of while loop

      
    %PH
    figure(1)
    pbaspect([1 1 1])
    semilogx(X_vec,PHss,'LineWidth',3)
    xlabel('Antagonist X)','FontSize',16,'FontWeight','bold')
    ylabel('PH','FontSize',16,'FontWeight','bold')
    xlim([1E-1 1E6])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northeast')
    legend boxon  

    %EH
    figure(2)
    pbaspect([1 1 1])
    semilogx(X_vec,EHss,'LineWidth',3)
    xlabel('Antagonist X','FontSize',16,'FontWeight','bold')
    ylabel('EH','FontSize',16,'FontWeight','bold')
    xlim([1E-1 1E6])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northeast')
    legend boxon  
    
    %EHCR
    figure(3)
    pbaspect([1 1 1])
    semilogx(X_vec,EHCRss,'LineWidth',3)
    xlabel('Antagonist X','FontSize',16,'FontWeight','bold')
    ylabel('EHCR','FontSize',16,'FontWeight','bold')
    xlim([1E-1 1E6])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northeast')
    legend boxon 

    %XCR
    figure(4)
    pbaspect([1 1 1])
    semilogx(X_vec,XCRss,'LineWidth',3)
    xlabel('Antagonist X','FontSize',16,'FontWeight','bold')
    ylabel('XCR','FontSize',16,'FontWeight','bold')
    xlim([1E-1 1E6])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northwest')
    legend boxon

    %EHPR
    figure(5)
    pbaspect([1 1 1])
    semilogx(X_vec,EHPRss,'LineWidth',3)
    xlabel('Antagonist X')
    ylabel('EHPR')
    xlim([1E-1 1E6])
    xticks([1E-2 1E0 1E2 1E4 1E6])
    set(gca,'FontSize',22)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northeast','FontSize',17)
    legend boxon

    %XPR
    figure(6)
    pbaspect([1 1 1])
    semilogx(X_vec,XPRss,'LineWidth',3)
    xlabel('Antagonist X')
    ylabel('XPR')
    xlim([1E-1 1E6])
    xticks([1E-2 1E0 1E2 1E4 1E6])
    set(gca,'FontSize',22)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northwest','FontSize',17)
    legend boxon
    
    %Endocrine Effect
    figure(7)
    pbaspect([2.2 1 1])
    semilogx(X_vec,Endocrine_Effect,'LineWidth',3)
    xlabel('Antagonist X')
    ylabel('Endocrine Effect')
    xlim([1E-1 1E6])
    xticks([1E-2 1E-1 1E0 1E1 1E2 1E3 1E4 1E5 1E6])
    set(gca,'FontSize',15)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northwest','FontSize',15)
    legend boxon 

    %Normalized Endocrine Effect
    figure(8)
    pbaspect([1 1 1])
    semilogx(X_vec,Normalized_Endocrine_Effect,'LineWidth',3)
    xlabel('Antagonist X','FontSize',16,'FontWeight','bold')
    ylabel('Normalized Endocrine Effect','FontSize',16,'FontWeight','bold')
    xlim([1E-1 1E6])
    xticks([1E-2 1E-1 1E0 1E1 1E2 1E3 1E4 1E5 1E6])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northwest')
    legend boxon 
    
    %EE, EH, and PR
    figure(9)
    pbaspect([1.5 1 1])
    loglog(X_vec,Endocrine_Effect,'LineWidth',4,'color',bluetints(:,j))%[0.00,0.45,0.74])
    hold on
    %A different legend name, legendInfo3, is used to avoid mix-up with the
    %name above because here we are plotting 3 lines at a time.
    legendInfo3{(j-1)*3+1} = [char(param_name), ' x ', num2str(fold_change_vec(j))];

    loglog(X_vec,EHss,':','LineWidth',3,'color',redtints(:,j))
    legendInfo3{(j-1)*3+2} = [char(param_name), ' x ', num2str(fold_change_vec(j))];

    loglog(X_vec,PRss,':','LineWidth',3,'color',greentints(:,j))
    legendInfo3{(j-1)*3+3} = [char(param_name), ' x ', num2str(fold_change_vec(j))];

    xlabel('Antagonist X','FontSize',16)
    ylabel('EE, EH, and PR','FontSize',16)
    xlim([1E-1 1E6])
    xticks([1E-2 1E-1 1E0 1E1 1E2 1E3 1E4 1E5 1E6])
    set(gca,'FontSize',12)
    box off
    legend(legendInfo3,'location','northeast')
    legend boxon

    %CR
    figure(10)
    pbaspect([1 1 1])
    semilogx(X_vec,CRss,'LineWidth',3)
    xlabel('Antagonist X','FontSize',16,'FontWeight','bold')
    ylabel('CR','FontSize',16,'FontWeight','bold')
    xlim([1E-1 1E6])
    xticks([1E-2 1E0 1E2 1E4 1E6])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northwest')
    legend boxon
    
    %PR
    figure(11)
    pbaspect([1 1 1])
    semilogx(X_vec,PRss,'LineWidth',3)
    xlabel('Antagonist X','FontSize',16,'FontWeight','bold')
    ylabel('PR','FontSize',16,'FontWeight','bold')
    xlim([1E-1 1E6])
    xticks([1E-2 1E0 1E2 1E4 1E6])
    set(gca,'FontSize',12)
    box off
    legendInfo{j} = [char(param_name), ' x ', num2str(fold_change_vec(j))];
    hold on
    legend(legendInfo,'location','northwest')
    legend boxon

    
end %End of j for loop



%% ------------------------- 6-PARAMETER MC SIMULATION -------------------------------------------------------------------------------%%

num_of_runs = 20000;

kf5_vec = [];
kf6_vec = [];
kf7_vec = [];
kf8_vec = [];
wp_vec = [];
wc_vec = [];

X_vec = [];
EE_curve = [];


for j = 1:1:num_of_runs

    j
    
    %-----------------Randomize the 6 parameters---------------%
    param.kf5 = 10^(-1+2*rand()) * default_param.kf5; %10^(-1+2*rand()) generates uniform distribution between 0.1-10 on log10 scale
    param.kf6 = 10^(-1+2*rand()) * default_param.kf6; 
    param.kf7 = 10^(-1+2*rand()) * default_param.kf7; 
    param.kf8 = 10^(-1+2*rand()) * default_param.kf8;
    
    param.wp = rand();
    param.wc = rand(); 

    
    %------------------Store parameters------------------%
    kf5_vec(j) = param.kf5;
    kf6_vec(j) = param.kf6;
    kf7_vec(j) = param.kf7;
    kf8_vec(j) = param.kf8;
    wp_vec(j)  = param.wp;
    wc_vec(j)  = param.wc;
    
       
    %EDC concentrations
    param.X = 0;
    
    i = 1;
    while (param.X <= X_final)  

        options2 = odeset('RelTol',1e-8,'AbsTol',1e-20);
        [t,y] = ode23tb('EDC_ode',tspan, y0, options2, param); 

        %Steay-state values of model output
        EHss(i)=y(end,1);
        PHss(i)=y(end,2);
        EHCRss(i)=y(end,3);
        XCRss(i)=y(end,4);
        EHPRss(i)=y(end,5);
        XPRss(i)=y(end,6);
        CRss(i)=param.CRtot-EHCRss(i)-XCRss(i);
        PRss(i)=param.PRtot-EHPRss(i)-XPRss(i);
        Endocrine_Effect(i)=EHPRss(i)+param.wp*XPRss(i);
        Endocrine_Effect_BasalValue_vec(j) = Endocrine_Effect(1);
        Normalized_Endocrine_Effect(i) = Endocrine_Effect(i)/Endocrine_Effect(1);
        X_vec(i) = param.X;
        

        %Increase X concentration
        if (param.X == 0)
            param.X = X_ini;
        else
            param.X = param.X * X_interval;
        end
        
        i = i + 1;

    end %End of while loop

    %Store the EE curves
    EE_curve(j,:) = Endocrine_Effect;
    
end %End of j loop

save('6_parameter_MC_simulation_results/MC_simulation_results_6_params_w_0-1_Uniform.mat', 'X_vec', 'EE_curve', 'kf5_vec','kf6_vec', 'kf7_vec', 'kf8_vec', 'wp_vec', 'wc_vec')

%% ------------------------- 6-PARAMETER MC SIMULATION RESULT ANALYSIS (Figs. 7-8, Figs. S8, S10A-S10D) ------------------------------%%

%Load previously saved DR curves
load('6_parameter_MC_simulation_results/MC_simulation_results_6_params_w_0-1_Uniform.mat');


% ---------- Group EE curves into different Monotonic and NMDR curves--------------%
EE_curve_Flat = [];
EE_curve_Increasing = [];
EE_curve_Decreasing = [];
EE_curve_U = [];
EE_curve_Bell = [];
EE_curve_U_then_Bell = [];
EE_curve_Bell_then_U = [];

num_of_curves_to_display = 50;

num_curves = length(EE_curve(:,1)); %number of curves
for i = 1:1:num_curves
    % calculate the difference between two neighboring points on each EE curve
    EE_curve_difference = diff(EE_curve(i,:));
    
    up_segment = any(EE_curve_difference>0);    %Boolean value, evaluate if the curve contains any positive slope
    down_segment = any(EE_curve_difference<0);  %Boolean value, evaluate if the curve contains any negative slope
    
    %Flat curve
    if not(up_segment) && not(down_segment) 
        EE_curve_Flat = [EE_curve_Flat; [i,EE_curve(i,:)]]; %Save the index of the curve and the curve itself
    
    %Monotonic increasing
    elseif up_segment && not(down_segment)   
        EE_curve_Increasing = [EE_curve_Increasing; [i,EE_curve(i,:)]];
        
    %Monotonic decreasing    
    elseif not(up_segment) && down_segment   
        EE_curve_Decreasing = [EE_curve_Decreasing; [i,EE_curve(i,:)]];
    
    %NMDR    
    else %up_segment && down_segment
        slope_sign_changes = diff(sign(EE_curve_difference)); %slope sign change between two neighboring segments
        
        %U shape
        if min(slope_sign_changes) >= 0
            EE_curve_U = [EE_curve_U; [i,EE_curve(i,:)]];
            
        %Bell shape
        elseif max(slope_sign_changes) <= 0
            EE_curve_Bell = [EE_curve_Bell; [i,EE_curve(i,:)]];
            
        %U_then_Bell shape. May contain U-Bell-U....
        elseif find(slope_sign_changes>0,1) < find(slope_sign_changes<0,1) %when the first positive sign change is before the first negative sign change
            EE_curve_U_then_Bell = [EE_curve_U_then_Bell; [i,EE_curve(i,:)]];
        
        %Bell_then_U shape. May contain Bell-U-Bell....
        else %find(slope_sign_changes>0,1) > find(slope_sign_changes<0,1) %when the first positive sign change is after the first negative sign change
            EE_curve_Bell_then_U = [EE_curve_Bell_then_U; [i,EE_curve(i,:)]];  
        end
                 
    end
end
    


% Plot flat curves
if not(isempty(EE_curve_Flat))
    figure(2100)
    semilogx(X_vec,EE_curve_Flat(1:min(num_of_curves_to_display,size(EE_curve_Flat,1)), 2:end),'LineWidth',2) %The first column is the index in EE_curve
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("Flat")
    box off
end

% Fig. S8A: Plot monotonically increasing curves
if not(isempty(EE_curve_Increasing))
    figure(2200)
    semilogx(X_vec,EE_curve_Increasing(1:min(num_of_curves_to_display,size(EE_curve_Increasing,1)), 2:end),'LineWidth',2)
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("MI")
    box off
end

% Fig. S8C: Plot monotonically decreasing curves
if not(isempty(EE_curve_Decreasing))
    figure(2300)
    semilogx(X_vec,EE_curve_Decreasing(1:min(num_of_curves_to_display,size(EE_curve_Decreasing,1)), 2:end),'LineWidth',2)
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("MD")
    box off
end

% Fig. S8B: Plot U-shape curves
if not(isempty(EE_curve_U))
    figure(2400)
    semilogx(X_vec,EE_curve_U(1:min(num_of_curves_to_display,size(EE_curve_U,1)) ,2:end),'LineWidth',2)
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("J/U")
    box off
end

% Fig. S8D: Plot Bell-shape curves
if not(isempty(EE_curve_Bell))
    figure(2500)
    semilogx(X_vec,EE_curve_Bell(1:min(num_of_curves_to_display,size(EE_curve_Bell,1)) ,2:end),'LineWidth',2)
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("Bell")
    box off
end

% Plot U-then-Bell-shape curves
if not(isempty(EE_curve_U_then_Bell))
    figure(2600)
    semilogx(X_vec,EE_curve_U_then_Bell(1:min(num_of_curves_to_display,size(EE_curve_U_then_Bell,1)) ,2:end),'LineWidth',2)
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("U-then-Bell")
    box off
end

% Plot Bell-then-U-shape curves
if not(isempty(EE_curve_Bell_then_U))
    figure(2700)
    semilogx(X_vec,EE_curve_Bell_then_U(1:min(num_of_curves_to_display,size(EE_curve_Bell_then_U,1)) ,2:end),'LineWidth',3)
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("Bell-then-U")
    box off
end


% Fig. S8E: Pie chart 
figure(2800)
Curve_numbers = [size(EE_curve_Increasing, 1), size(EE_curve_U, 1), size(EE_curve_Decreasing, 1),  size(EE_curve_Bell, 1), size(EE_curve_U_then_Bell, 1), size(EE_curve_Bell_then_U,1)];
h = pie(Curve_numbers)
set(findobj(h,'type','text'),'fontsize',18)
labels = {'MI','J/U','MD','Bell'};
legend(labels, 'Location','eastoutside', 'FontSize',16);



% ------------Analyze parameter combinations that can differentiate different curve shapes-----%
II = EE_curve_Increasing(:,1);
UI = EE_curve_U(:,1);
DI = EE_curve_Decreasing(:,1);
BI = EE_curve_Bell(:,1);

Kd5_II = default_param.kb5 ./ kf5_vec(II);
Kd6_II = default_param.kb6 ./ kf6_vec(II);
Kd7_II = default_param.kb7 ./ kf7_vec(II);
Kd8_II = default_param.kb8 ./ kf8_vec(II);

Kd5_UI = default_param.kb5 ./ kf5_vec(UI);
Kd6_UI = default_param.kb6 ./ kf6_vec(UI);
Kd7_UI = default_param.kb7 ./ kf7_vec(UI);
Kd8_UI = default_param.kb8 ./ kf8_vec(UI);

Kd5_DI = default_param.kb5 ./ kf5_vec(DI);
Kd6_DI = default_param.kb6 ./ kf6_vec(DI);
Kd7_DI = default_param.kb7 ./ kf7_vec(DI);
Kd8_DI = default_param.kb8 ./ kf8_vec(DI);

Kd5_BI = default_param.kb5 ./ kf5_vec(BI);
Kd6_BI = default_param.kb6 ./ kf6_vec(BI);
Kd7_BI = default_param.kb7 ./ kf7_vec(BI);
Kd8_BI = default_param.kb8 ./ kf8_vec(BI);

wp_II = wp_vec(II);
wc_II = wc_vec(II);
wp_UI = wp_vec(UI);
wc_UI = wc_vec(UI);

wp_DI = wp_vec(DI);
wc_DI = wc_vec(DI);
wp_BI = wp_vec(BI);
wc_BI = wc_vec(BI);


% Fig. 7A: Scatter plot for parameter conditions differentiating monotonic increasing vs. U-shape curves (1000 randomly selected)
x_II = Kd5_II ./ Kd6_II .* wp_II;
y_II = Kd7_II ./ Kd8_II .* wc_II;

x_UI = Kd5_UI ./ Kd6_UI .* wp_UI;
y_UI = Kd7_UI ./ Kd8_UI .* wc_UI;

figure(1001)
scatter(log10(x_II(1:1000)), log10(y_II(1:1000)), 10, 'filled')
legendInfo{1} = ['MI'];
hold on
scatter(log10(x_UI(1:1000)), log10(y_UI(1:1000)), 10, 'filled', 'MarkerFaceAlpha',0.8)
legendInfo{2} = ['J/U'];
xlim([-4,3])
xticks([-4:1:3])
ylim([-4,3])
yticks([-4:1:3])
xlabel("Log10 (Kd5/Kd6*Wp)")
ylabel("Log10 (Kd7/Kd8*Wc)")
pbaspect([1 1 1])
set(gca,'FontSize',18)
box off
x = [-4:3];
y = [-4:3];
plot(x,y, 'black', 'HandleVisibility', 'off')
legend(legendInfo,'location','northwest')
legend boxon  
title("MI vs. J/U")


% Fig. 7B: Histograms for parameter conditions differentiating monotonic increasing vs. U-shape curves
composite_II = (Kd7_II ./ Kd8_II .* wc_II) ./ (Kd5_II ./ Kd6_II .* wp_II);
composite_UI = (Kd7_UI ./ Kd8_UI .* wc_UI) ./ (Kd5_UI ./ Kd6_UI .* wp_UI);
figure (1002)
hold on
histogram(log10(composite_II), 'BinWidth', 0.5, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(composite_UI), 'BinWidth', 0.5, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(composite_II), 'BinWidth', 0.5, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
histogram(log10(composite_UI), 'BinWidth', 0.5, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
xlim([-5,5])
xticks([-4:1:4])
xlabel("Log10 ((Kd7/Kd8*Wc) / (Kd5/Kd6*Wp))")
ylabel('Frequency')
pbaspect([1 1 1])
set(gca,'FontSize',18)
legend(legendInfo,'location','northwest')
legend boxon  
title("Increasing vs. J/U-shape")


% Fig. 7C: U-shape curve magnitude vs Log10 ((Kd7/Kd8*Wc) / (Kd5/Kd6*Wp))
MagL = EE_curve_U(:,2) - min(EE_curve_U(:,2:end),[], 2);
MagL_normalized = (EE_curve_U(:,2) - min(EE_curve_U(:,2:end),[], 2)) ./ EE_curve_U(:,2);

MagR = EE_curve_U(:,end) - min(EE_curve_U(:,2:end),[], 2);
MagR_normalized = (EE_curve_U(:,end) - min(EE_curve_U(:,2:end),[], 2)) ./ EE_curve_U(:,end);

figure(1003)
scatter(log10(composite_UI), MagL_normalized, 10, 'filled')
xlabel("Log10 ((Kd7/Kd8*Wc) / (Kd5/Kd6*Wp))")
ylabel('Magnitude of J/U')
xlim([-2,4])
xticks([-2:1:4])
pbaspect([1 1 1])
set(gca,'FontSize',18)


% Fig. 7D: Scatter plot for parameter conditions differentiating monotonic decreasing vs. Bell-shape curves
figure(1004)
scatter(log10(Kd8_DI ./ Kd5_DI .* wc_DI), log10(Kd6_DI ./ Kd7_DI .* wp_DI), 10, 'filled')
legendInfo{1} = ['MD'];
hold on
scatter(log10(Kd8_BI ./ Kd5_BI .* wc_BI), log10(Kd6_BI ./ Kd7_BI .* wp_BI), 10, 'filled', 'MarkerFaceAlpha',0.8)
legendInfo{2} = ['Bell'];
xlim([-5,3])
xticks([-5:1:3])
ylim([-5,3])
yticks([-5:1:3])
xlabel("Log10 (Kd8/Kd5*Wc)")
ylabel("Log10 (Kd6/Kd7*Wp)")
pbaspect([1 1 1])
set(gca,'FontSize',18)
box off
x = [-5:3];
y = [-5:3];
plot(x,y, 'black', 'HandleVisibility', 'off')
legend(legendInfo,'location','northwest')
legend boxon  
title("MD vs. Bell")


% Fig. 7E: Histograms for parameter conditions differentiating monotonic decreasing vs. Bell-shape curves
composite_DI = (Kd6_DI ./ Kd7_DI .* wp_DI) ./ (Kd8_DI ./ Kd5_DI .* wc_DI);
composite_BI = (Kd6_BI ./ Kd7_BI .* wp_BI) ./ (Kd8_BI ./ Kd5_BI .* wc_BI);
figure (1005)
hold on
histogram(log10(composite_DI), 'BinWidth', 0.5, 'EdgeColor', 'none', 'DisplayName','MD')
histogram(log10(composite_BI), 'BinWidth', 0.5, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
histogram(log10(composite_DI), 'BinWidth', 0.5, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
histogram(log10(composite_BI), 'BinWidth', 0.5, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
xlim([-6,5])
xticks([-6:2:4])
xlabel("Log10 ((Kd6/Kd7*Wp) / (Kd8/Kd5*Wc))")
ylabel('Frequency')
pbaspect([1 1 1])
set(gca,'FontSize',18)
legend(legendInfo,'location','northwest')
legend boxon
title("Decreasing vs. Bell-shape")


% Fig. 7F: Bell-shape curve magnitude vs Log10 ((Kd6/Kd7*Wp) / (Kd8/Kd5*Wc))
MagL = max(EE_curve_Bell(:,2:end),[], 2) - EE_curve_Bell(:,2);
MagL_normalized = (max(EE_curve_Bell(:,2:end),[], 2) - EE_curve_Bell(:,2)) ./ EE_curve_Bell(:,2);

MagR = max(EE_curve_Bell(:,2:end),[], 2) - EE_curve_Bell(:,end);
MagR_normalized = (max(EE_curve_Bell(:,2:end),[], 2) - EE_curve_Bell(:,end)) ./ EE_curve_Bell(:,end);

figure(1006)
scatter(log10(composite_BI), log10(MagL_normalized), 10, 'filled')
xlabel("Log10 ((Kd6/Kd7*Wp) / (Kd8/Kd5*Wc))")
ylabel('Magnitude of Bell')
xlim([-2,4])
xticks([-2:1:4])
pbaspect([1 1 1])
set(gca,'FontSize',18)



% ------------Histograms for Wp and Wc for differentiating MI vs MD curves---------------%
% Fig. S10A: Histograms for Wp for differentiating MI vs MD curves
figure (10041)
hold on
histogram((wp_II), 'BinWidth', 0.02, 'EdgeColor', 'none', 'DisplayName','MI')
histogram((wp_DI), 'BinWidth', 0.02, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram((wp_II), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
histogram((wp_DI), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
xticks([0, 0.1, 0.5, 1])
xlabel("Wp")
ylabel('Frequency')
pbaspect([1 1 1])
set(gca,'FontSize',18)
pbaspect([1.4 1 1])
legend('location','northwest')
legend boxon  
title("MI vs. MD")

% Fig. S10B: Histograms for Wc for differentiating MI vs MD curves
figure (10042)
hold on
histogram((wc_II), 'BinWidth', 0.02, 'EdgeColor', 'none', 'DisplayName','MI')
histogram((wc_DI), 'BinWidth', 0.02, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram((wc_II), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
histogram((wc_DI), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
xticks([0, 0.1, 0.5, 1])
xlabel("Wc")
ylabel('Frequency')
pbaspect([1 1 1])
set(gca,'FontSize',18)
pbaspect([1.4 1 1])
legend('location','northwest')
legend boxon  
title("MI vs. MD")


% ------------Histograms for Wc and Wp for differentiating U-shape vs Bell-shape curves---------------%
% Fig. S10C: Histograms for Wp for differentiating U-shape vs Bell-shape curves
figure (10043)
hold on
histogram((wp_UI), 'BinWidth', 0.02, 'EdgeColor', 'none', 'DisplayName','J/U')
histogram((wp_BI), 'BinWidth', 0.02, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
histogram((wp_UI), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
histogram((wp_BI), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
xticks([0, 0.1, 0.5, 1])
xlabel("Wp")
ylabel('Frequency')
pbaspect([1 1 1])
set(gca,'FontSize',18)
pbaspect([1.4 1 1])
legend('location','northwest')
legend boxon  
title("J/U-shape vs. Bell-shape")

% Fig. S10D: Histograms for Wc for differentiating U-shape vs Bell-shape curves
figure (10044)
hold on
histogram((wc_UI), 'BinWidth', 0.02, 'EdgeColor', 'none', 'DisplayName','J/U')
histogram((wc_BI), 'BinWidth', 0.02, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
histogram((wc_UI), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
histogram((wc_BI), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
xticks([0, 0.1, 0.5, 1])
xlabel("Wc")
ylabel('Frequency')
pbaspect([1 1 1])
set(gca,'FontSize',18)
pbaspect([1.4 1 1])
legend('location','northwest')
legend boxon  
title("J/U-shape vs. Bell-shape")


% % Histograms for wp/wc for differentiating U-shape vs Bell-shape curves
% composite_UI = wc_UI ./ wp_UI;
% composite_BI = wc_BI ./ wp_BI;
% figure (10045)
% hold on
% histogram(log10(composite_UI), 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','J/U')
% histogram(log10(composite_BI), 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
% histogram(log10(composite_UI), 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
% histogram(log10(composite_BI), 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
% xlabel("Log10 (Wc / Wp)")
% ylabel('Frequency')
% pbaspect([1 1 1])
% set(gca,'FontSize',18)
% legend('location','northwest')
% legend boxon  
% title("J/U-shape vs. bell-shape")
% 
% 
% % Histograms for Kd7*wc for differentiating U-shape vs bell-shape curves
% composite_UI = Kd7_UI .* wc_UI;
% composite_BI = Kd7_BI .* wc_BI;
% figure (10046)
% hold on
% histogram(log10(composite_UI), 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','J/U')
% histogram(log10(composite_BI), 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
% histogram(log10(composite_UI), 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
% histogram(log10(composite_BI), 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
% %xlim([-5,5])
% %xticks([-4:1:4])
% xlabel("Log10 (Kd7 * Wc)")
% ylabel('Frequency')
% pbaspect([1 1 1])
% set(gca,'FontSize',18)
% legend('location','northwest')
% legend boxon  
% title("J/U-shape vs. bell-shape")


% %-------------------------Scatter plot for individual parameter pairs among the 4 curve types------------------%
% figure(10047)
% scatter(log10(Kd7_DI), log10(Kd6_DI), 10, 'filled')
% legendInfo{1} = ['MD'];
% hold on
% scatter(log10(Kd7_BI), log10(Kd6_BI), 10, 'filled', 'MarkerFaceAlpha',0.8)
% legendInfo{2} = ['Bell'];
% scatter(log10(Kd7_II), log10(Kd6_II), 10, 'filled', 'MarkerFaceAlpha',0.8)
% legendInfo{3} = ['MI'];
% scatter(log10(Kd7_UI), log10(Kd6_UI), 10, 'filled', 'MarkerFaceAlpha',0.8)
% legendInfo{4} = ['UI'];
% % xlim([-5,3])
% % xticks([-5:1:3])
% % ylim([-5,3])
% % yticks([-5:1:3])
% xlabel("Log10 (Kd7)")
% ylabel("Log10 (Kd6)")
% pbaspect([1 1 1])
% set(gca,'FontSize',18)
% box off
% % x = [-5:3];
% % y = [-5:3];
% % plot(x,y, 'black', 'HandleVisibility', 'off')
% legend(legendInfo,'location','northwest')
% legend boxon  
% title("Kd7 vs. Kd6")


%----------------Histogram of individual parameters for the 4 curve shapes----------------------%

% Fig. 8A
default_Kd5 = default_param.kb5/default_param.kf5;
% Kd5 - MI and J/U
figure(2005)
hold on
histogram(log10(Kd5_II/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(Kd5_UI/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(Kd5_II/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd5_UI/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (Kd5 fold change)")
ylabel('Frequency (x1000)')
ylim([0,1000])
yticks([0:500:1000])
yticklabels({'0', '0.5', '1'})
pbaspect([1 1 1])
set(gca,'FontSize',24)
legend boxoff
legend off

% Kd5 - MD and Bell
figure(20051)
hold on
histogram(log10(Kd5_DI/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd5_BI/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd5_DI/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(Kd5_BI/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (Kd5 fold change)")
ylabel('Frequency (x1000)')
yticks([0:50:100])
yticklabels({'0', '0.05', '0.1'})
pbaspect([1 1 1])
set(gca,'FontSize',24)
legend boxoff
legend off


% Fig. 8B
default_Kd6 = default_param.kb6/default_param.kf6;
% Kd6 - MI and J/U
figure(2006)
hold on
histogram(log10(Kd6_II/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(Kd6_UI/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(Kd6_II/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd6_UI/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (Kd6 fold change)")
ylabel('Frequency (x1000)')
ylim([0,1000])
yticks([0:500:1000])
yticklabels({'0', '0.5', '1'})
pbaspect([1 1 1])
set(gca,'FontSize',24)
legend boxoff
legend off


% Kd6 - MD and Bell
figure(20061)
hold on
histogram(log10(Kd6_DI/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd6_BI/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd6_DI/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(Kd6_BI/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (Kd6 fold change)")
ylabel('Frequency (x1000)')
ylim([0,100])
yticks([0:50:100])
yticklabels({'0', '0.05', '0.1'})
pbaspect([1 1 1])
set(gca,'FontSize',24)
legend boxoff  
legend off


% Fig. 8C
default_Kd7 = default_param.kb7/default_param.kf7;
% Kd7 - MI and J/U
figure(2007)
hold on
histogram(log10(Kd7_II/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(Kd7_UI/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(Kd7_II/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd7_UI/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (Kd7 fold change)")
ylabel('Frequency (x1000)')
ylim([0,1000])
yticks([0:500:1000])
yticklabels({'0', '0.5', '1'})
pbaspect([1 1 1])
set(gca,'FontSize',24)
legend boxoff
legend off

% Kd7 - MD and Bell
figure(20071)
hold on
histogram(log10(Kd7_DI/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd7_BI/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd7_DI/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(Kd7_BI/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (Kd7 fold change)")
ylabel('Frequency (x1000)')
ylim([0,100])
yticks([0:50:100])
yticklabels({'0', '0.05', '0.1'})
pbaspect([1 1 1])
set(gca,'FontSize',24)
legend boxoff  
legend off


% Fig. 8D
default_Kd8 = default_param.kb8/default_param.kf8;
% Kd8 - MI and J/U
figure(2008)
hold on
histogram(log10(Kd8_II/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(Kd8_UI/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(Kd8_II/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd8_UI/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (Kd8 fold change)")
ylabel('Frequency (x1000)')
ylim([0,1000])
yticks([0:500:1000])
yticklabels({'0', '0.5', '1'})
pbaspect([1 1 1])
set(gca,'FontSize',24)
legend boxoff
legend off

% Kd8 - MD and Bell
figure(20081)
hold on
histogram(log10(Kd8_DI/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd8_BI/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd8_DI/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(Kd8_BI/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (Kd8 fold change)")
ylabel('Frequency (x1000)')
ylim([0,100])
yticks([0:50:100])
yticklabels({'0', '0.05', '0.1'})
pbaspect([1 1 1])
set(gca,'FontSize',24)
legend boxoff  
legend off


% Fig. 8E
default_wp = default_param.wp;
% wp - MI and J/U
figure(2009)
hold on
histogram(wp_II/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(wp_UI/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(wp_II/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(wp_UI/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("wp")
ylabel('Frequency (x1000)')
ylim([0,1000])
yticks([0:500:1000])
yticklabels({'0', '0.5', '1'})
pbaspect([1 1 1])
set(gca,'FontSize',24)
legend boxoff
legend off

% wp - MD and Bell
figure(20091)
hold on
histogram(wp_DI/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(wp_BI/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(wp_DI/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(wp_BI/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("wp")
ylabel('Frequency (x1000)')
ylim([0,400])
yticks([0:100:400])
yticklabels({'0', '0.1', '0.2','0.3','0.4'})
pbaspect([1 1 1])
set(gca,'FontSize',24)
legend boxoff  
legend off


% Fig. 8F
default_wc = default_param.wc;
% wc - MI and J/U
figure(2010)
hold on
histogram(wc_II/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(wc_UI/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(wc_II/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(wc_UI/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("wc")
ylabel('Frequency (x1000)')
ylim([0,1000])
yticks([0:500:1000])
yticklabels({'0', '0.5', '1'})
pbaspect([1 1 1])
set(gca,'FontSize',24)
legend boxoff
legend off


% wc - MD and Bell
figure(20101)
hold on
histogram(wc_DI/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(wc_BI/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(wc_DI/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(wc_BI/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("wc")
ylabel('Frequency (x1000)')
ylim([0,500])
yticks([0:100:500])
yticklabels({'0', '0.1', '0.2','0.3','0.4','0.5'})
pbaspect([1 1 1])
set(gca,'FontSize',24)
legend boxoff  
legend off

%% ------------------------- POPULATION MC SIMULATION---------------------------------------------------------------------------------%%

% Load population parameters;
%load('Population_MC_simulation_results/parameter_matrix_10xk3_300xk1.mat');
load('Population_MC_simulation_results/100xk3_1000xk1/param_9996.mat');
parameter_matrix = param_matrix';

% Load population steady-state state variables;
%load('Population_MC_simulation_results/300xk3_1000xk1/hormone_9996.mat');

EE_curve = [];
EH_curve = [];

num_of_runs = length(parameter_matrix);

for j = 1:1:num_of_runs

    j
      
    %-----------------Obtain parameters from the population parameter matrix---------------%
    param.n3        = parameter_matrix(j,1);
    param.kf7       = parameter_matrix(j,2);
    param.k30       = parameter_matrix(j,3);
    param.k3        = parameter_matrix(j,4);
    param.k1        = parameter_matrix(j,5);
    param.Kd3       = parameter_matrix(j,6);
    param.CRtot     = parameter_matrix(j,7);
    
    %------------------Randomize relevant receptor binding and efficacy parameters------------------%
    param.kf5 = 10^(-1+2*rand())*default_param.kf5; %10^(-1+2*rand()) generates uniform distribution between 0.1-10 on log10 scale
    param.kf6 = 10^(-1+2*rand())*default_param.kf6; 
    param.kf8 = 10^(-1+2*rand())*default_param.kf8;
    
    param.wp = rand();
    param.wc = rand(); 

    %PRtot is not varied since it will have no effect on the result.
    
%    %generate wp and wc that have equal chance to be between 0-0.1 and 0.1-1.0
%     p = rand();
%     if p <0.5
%         param.wp = 0.1*rand()*default_param.wp;
%     else
%         param.wp = (0.1+0.9*rand())*default_param.wp; %uniform distribution between 0.1-1.0
%     end
%     
%     c = rand();
%     if c <0.5
%         param.wc = 0.1*rand()*default_param.wc;
%     else
%         param.wc = (0.1+0.9*rand())*default_param.wc; %uniform distribution between 0.1-1.0
%     end    
    
    
    %------------------Store parameters------------------%
    kf5_vec(j) = param.kf5;
    kf6_vec(j) = param.kf6;
    kf7_vec(j) = param.kf7;
    kf8_vec(j) = param.kf8;
    wp_vec(j) = param.wp;
    wc_vec(j) = param.wc;

    n3_vec(j) = param.n3;
    k30_vec(j) = param.k30;
    k3_vec(j) = param.k3;
    k1_vec(j) = param.k1;
    Kd3_vec(j) = param.Kd3;
    CRtot_vec(j) = param.CRtot;
    
        
    %EDC concentrations
    param.X = 0;
    

    i = 1;
    while (param.X <= X_final)  

        %make sure the solution will not fall within the negative range
        options2 = odeset('RelTol',1e-8,'AbsTol',1e-20);
        %after multiple tests, ode23tb solver is more efficient compared to other solvers 
        [t,y] = ode23tb('EDC_ode',tspan, y0, options2, param); 

        %Steay-state values of model output
        EHss(i)=y(end,1);
        PHss(i)=y(end,2);
        EHCRss(i)=y(end,3);
        XCRss(i)=y(end,4);
        EHPRss(i)=y(end,5);
        XPRss(i)=y(end,6);
        CRss(i)=param.CRtot-EHCRss(i)-XCRss(i);
        PRss(i)=param.PRtot-EHPRss(i)-XPRss(i);
        Endocrine_Effect(i)=EHPRss(i)+param.wp*XPRss(i);
        Endocrine_Effect_BasalValue_vec(j) = Endocrine_Effect(1);
        Normalized_Endocrine_Effect(i) = Endocrine_Effect(i)/Endocrine_Effect(1);
        X_vec(i) = param.X;
        

        %Increase X concentration
        if (param.X == 0)
            param.X = X_ini;
        else
            param.X = param.X * X_interval;
        end
        
        i = i + 1;

    end %End of while loop

    
    
    
    %Store the EE curves and others
    EE_curve(j,:) = Endocrine_Effect;
    EH_curve(j,:) = EHss;

end %End of j loop

save('Population_MC_simulation_results/100xk3_1000xk1/Population_MC_simulation_results.mat', 'X_vec', 'EE_curve', 'EH_curve', 'kf5_vec','kf6_vec', 'kf7_vec', 'kf8_vec', 'wp_vec', 'wc_vec', 'n3_vec', 'k30_vec', 'k3_vec', 'k1_vec', 'Kd3_vec', 'CRtot_vec')

toc

%% ------------------------- POPULATION MC SIMULATION RESULT ANALYSIS (Figs. 9-10, Figs. S9, S10E-S10H) ------------------------------%%

load('Population_MC_simulation_results/hormone_9996.mat');
hormone_matrix = hormone_matrix';

load('Population_MC_simulation_results/param_9996.mat');
parameter_matrix = param_matrix';

%%Load previously saved DR curves
load('Population_MC_simulation_results/Population_MC_simulation_results.mat');


% %---------------Display baseline population hormone and parameter histograms------------%
% %Scatter plot of PH and EH
% figure(1000)
% scatter(hormone_matrix(:,1), hormone_matrix(:,2), '.')
% set(gca, 'YScale', 'log')
% xlabel('EH')
% ylabel('log10(PH)')
% pbaspect([1 1 1])
% set(gca,'FontSize',18)
% box off
% 
% %EH
% mean_EH     = mean(hormone_matrix(:,1));
% std_EH      = std(hormone_matrix(:,1));
% figure(1001)
% histogram(hormone_matrix(:,1))
% xlabel('EH')
% title(strcat('mean=', num2str(mean_EH), ', std=', num2str(std_EH)))
% 
% %PH
% mean_PH     = mean(hormone_matrix(:,2));
% std_PH      = std(hormone_matrix(:,2));
% figure(1002)
% histogram(log10(hormone_matrix(:,2)))
% xlabel('log10(PH)')
% title(strcat('mean=', num2str(mean_PH), ', std=', num2str(std_PH)))
% 
% %n3
% mean_n3     = mean(parameter_matrix(:,1));
% std_n3      = std(parameter_matrix(:,1));
% figure(1003)
% histogram(parameter_matrix(:,1))
% xlabel('n3')
% title(strcat('mean=', num2str(mean_n3), ', std=', num2str(std_n3)))
% hold on
% 
% %kf7
% mean_kf7     = mean(parameter_matrix(:,2));
% std_kf7      = std(parameter_matrix(:,2));
% figure(1004)
% histogram(log10(parameter_matrix(:,2)))
% xlabel('Log10 kf7')
% title(strcat('mean=', num2str(mean_kf7), ', std=', num2str(std_kf7)))
% hold on
% 
% %k30
% mean_k30     = mean(parameter_matrix(:,3));
% std_k30      = std(parameter_matrix(:,3));
% figure(1005)
% histogram(log10(parameter_matrix(:,3)))
% xlabel('Log10 k30')
% title(strcat('mean=', num2str(mean_k30), ', std=', num2str(std_k30)))
% hold on
% 
% %k3
% geomean_k3     = geomean(parameter_matrix(:,4));
% std_k3      = std(parameter_matrix(:,4));
% figure(1006)
% histogram(log10(parameter_matrix(:,4)))
% xlabel('Log10 k3')
% title(strcat('geomean=', num2str(geomean_k3), ', std=', num2str(std_k3)))
% hold on
% 
% %k1
% mean_k1     = mean(parameter_matrix(:,5));
% std_k1      = std(parameter_matrix(:,5));
% figure(1007)
% histogram(log10(parameter_matrix(:,5)))
% xlabel('Log10 k1')
% title(strcat('mean=', num2str(mean_k1), ', std=', num2str(std_k1)))
% hold on
% 
% %Kd3
% mean_Kd3     = mean(parameter_matrix(:,6));
% std_Kd3      = std(parameter_matrix(:,6));
% figure(1008)
% histogram(log10(parameter_matrix(:,6)))
% xlabel('Log10 Kd3')
% title(strcat('mean=', num2str(mean_Kd3), ', std=', num2str(std_Kd3)))
% hold on
% 
% %CRtot
% mean_CRtot     = mean(parameter_matrix(:,7));
% std_CRtot      = std(parameter_matrix(:,7));
% figure(1009)
% histogram(log10(parameter_matrix(:,7)))
% xlabel('Log10 CRtot')
% title(strcat('mean=', num2str(mean_CRtot), ', std=', num2str(std_CRtot)))
% hold on


% ---------- Group EE curves into different monotonic and NMDR curves--------------%
EE_curve_Flat = [];
EE_curve_Increasing = [];
EE_curve_Decreasing = [];
EE_curve_U = [];
EE_curve_Bell = [];
EE_curve_U_then_Bell = [];
EE_curve_Bell_then_U = [];

num_of_curves_to_display = 50;

num_curves = length(EE_curve(:,1)); %number of curves
for i = 1:1:num_curves
    % calculate the difference between two neighboring points on each EE curve
    EE_curve_difference = diff(EE_curve(i,:));
    
    up_segment = any(EE_curve_difference>0);    %Boolean value, evaluate if the curve contains any positive slope
    down_segment = any(EE_curve_difference<0);  %Boolean value, evaluate if the curve contains any negative slope
    
    %Flat curve
    if not(up_segment) && not(down_segment) 
        EE_curve_Flat = [EE_curve_Flat; [i,EE_curve(i,:)]]; %Save the index of the curve and the curve itself
    
    %Monotonic increasing
    elseif up_segment && not(down_segment)   
        EE_curve_Increasing = [EE_curve_Increasing; [i,EE_curve(i,:)]];
        
    %Monotonic decreasing    
    elseif not(up_segment) && down_segment   
        EE_curve_Decreasing = [EE_curve_Decreasing; [i,EE_curve(i,:)]];
    
    %NMDR    
    else %up_segment && down_segment
        slope_sign_changes = diff(sign(EE_curve_difference)); %slope sign change between two neighboring segments
        
        %U shape
        if min(slope_sign_changes) >= 0
            EE_curve_U = [EE_curve_U; [i,EE_curve(i,:)]];
            
        %Bell shape
        elseif max(slope_sign_changes) <= 0
            EE_curve_Bell = [EE_curve_Bell; [i,EE_curve(i,:)]];
            
        %U_then_Bell shape. May contain U-Bell-U....
        elseif find(slope_sign_changes>0,1) < find(slope_sign_changes<0,1) %when the first positive sign change is before the first negative sign change
            EE_curve_U_then_Bell = [EE_curve_U_then_Bell; [i,EE_curve(i,:)]];
        
        %Bell_then_U shape. May contain Bell-U-Bell....
        else %find(slope_sign_changes>0,1) > find(slope_sign_changes<0,1) %when the first positive sign change is after the first negative sign change
            EE_curve_Bell_then_U = [EE_curve_Bell_then_U; [i,EE_curve(i,:)]];  
        end
                 
    end
end
    

% Plot flat curves
if not(isempty(EE_curve_Flat))
    figure(2100)
    semilogx(X_vec,EE_curve_Flat(1:min(num_of_curves_to_display,size(EE_curve_Flat,1)), 2:end),'LineWidth',2) %The first column is the index in EE_curve
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("Flat")
    box off
end

% Fig. S9A: Plot monotonically increasing curves
if not(isempty(EE_curve_Increasing))
    figure(2200)
    semilogx(X_vec,EE_curve_Increasing(1:min(num_of_curves_to_display,size(EE_curve_Increasing,1)), 2:end),'LineWidth',2)
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("MI")
    box off
end

% Fig. S9C: Plot monotonically decreasing curves
if not(isempty(EE_curve_Decreasing))
    figure(2300)
    semilogx(X_vec,EE_curve_Decreasing(1:min(num_of_curves_to_display,size(EE_curve_Decreasing,1)), 2:end),'LineWidth',2)
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("MD")
    box off
end

% Fig. S9B: Plot U-shape curves
if not(isempty(EE_curve_U))
    figure(2400)
    semilogx(X_vec,EE_curve_U(1:min(num_of_curves_to_display,size(EE_curve_U,1)) ,2:end),'LineWidth',2)
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("J/U")
    box off
end

% Fig. S9D: Plot Bell-shape curves
if not(isempty(EE_curve_Bell))
    figure(2500)
    semilogx(X_vec,EE_curve_Bell(1:min(num_of_curves_to_display,size(EE_curve_Bell,1)) ,2:end),'LineWidth',2)
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("Bell")
    box off
end

% Fig. S9E: Plot U-then-Bell-shape curves
if not(isempty(EE_curve_U_then_Bell))
    figure(2600)
    semilogx(X_vec,EE_curve_U_then_Bell(1:min(num_of_curves_to_display,size(EE_curve_U_then_Bell,1)) ,2:end),'LineWidth',2)
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("U-then-Bell")
    box off
end

% Fig. S9F: Plot Bell-then-U-shape curves
if not(isempty(EE_curve_Bell_then_U))
    figure(2700)
    semilogx(X_vec,EE_curve_Bell_then_U(1:min(num_of_curves_to_display,size(EE_curve_Bell_then_U,1)) ,2:end),'LineWidth',3)
    pbaspect([1.5 1 1])
    set(gca,'FontSize',18)
    xlim([1e-2, 1e6])
    xticks([1e-2, 1e0, 1e2, 1e4, 1e6])
    ylim([0, 10])
    yticks([0:2:10])
    xlabel('X','FontSize',16)
    ylabel('EE','FontSize',16)
    title("Bell-then-U")
    box off
end

% Fig. S9G: Pie chart
figure(2800)
Curve_numbers = [size(EE_curve_Increasing, 1), size(EE_curve_U, 1), size(EE_curve_Decreasing, 1), size(EE_curve_Bell, 1), size(EE_curve_U_then_Bell, 1), size(EE_curve_Bell_then_U,1)];
h = pie(Curve_numbers)
set(findobj(h,'type','text'),'fontsize',18)
labels = {'MI', 'J/U', 'MD','Bell', 'U then Bell', 'Bell then U'};
legend(labels, 'Location','eastoutside', 'FontSize',16);



% ------------Analyze parameter combinations that can differentiate different curve shapes-----%
II = EE_curve_Increasing(:,1);
UI = EE_curve_U(:,1);
DI = EE_curve_Decreasing(:,1);
BI = EE_curve_Bell(:,1);

Kd5_II = default_param.kb5 ./ kf5_vec(II);
Kd6_II = default_param.kb6 ./ kf6_vec(II);
Kd7_II = default_param.kb7 ./ kf7_vec(II);
Kd8_II = default_param.kb8 ./ kf8_vec(II);

Kd5_UI = default_param.kb5 ./ kf5_vec(UI);
Kd6_UI = default_param.kb6 ./ kf6_vec(UI);
Kd7_UI = default_param.kb7 ./ kf7_vec(UI);
Kd8_UI = default_param.kb8 ./ kf8_vec(UI);

Kd5_DI = default_param.kb5 ./ kf5_vec(DI);
Kd6_DI = default_param.kb6 ./ kf6_vec(DI);
Kd7_DI = default_param.kb7 ./ kf7_vec(DI);
Kd8_DI = default_param.kb8 ./ kf8_vec(DI);

Kd5_BI = default_param.kb5 ./ kf5_vec(BI);
Kd6_BI = default_param.kb6 ./ kf6_vec(BI);
Kd7_BI = default_param.kb7 ./ kf7_vec(BI);
Kd8_BI = default_param.kb8 ./ kf8_vec(BI);

wp_II = wp_vec(II);
wc_II = wc_vec(II);
wp_UI = wp_vec(UI);
wc_UI = wc_vec(UI);

wp_DI = wp_vec(DI);
wc_DI = wc_vec(DI);
wp_BI = wp_vec(BI);
wc_BI = wc_vec(BI);


%Temporary use of parameter_matrix for vec
n3_II = parameter_matrix(II,1);
n3_UI = parameter_matrix(UI,1);
n3_DI = parameter_matrix(DI,1);
n3_BI = parameter_matrix(BI,1);

k30_II = parameter_matrix(II,3);
k30_UI = parameter_matrix(UI,3);
k30_DI = parameter_matrix(DI,3);
k30_BI = parameter_matrix(BI,3);

k3_II = parameter_matrix(II,4);
k3_UI = parameter_matrix(UI,4);
k3_DI = parameter_matrix(DI,4);
k3_BI = parameter_matrix(BI,4);

k1_II = parameter_matrix(II,5);
k1_UI = parameter_matrix(UI,5);
k1_DI = parameter_matrix(DI,5);
k1_BI = parameter_matrix(BI,5);

Kd3_II = parameter_matrix(II,6);
Kd3_UI = parameter_matrix(UI,6);
Kd3_DI = parameter_matrix(DI,6);
Kd3_BI = parameter_matrix(BI,6);

CRtot_II = parameter_matrix(II,7);
CRtot_UI = parameter_matrix(UI,7);
CRtot_DI = parameter_matrix(DI,7);
CRtot_BI = parameter_matrix(BI,7);

PRtot_II = parameter_matrix(II,7);
PRtot_UI = parameter_matrix(UI,7);
PRtot_DI = parameter_matrix(DI,7);
PRtot_BI = parameter_matrix(BI,7);



% Fig. 9A: Scatter plot for parameter conditions differentiating monotonic increasing vs. U-shape curves (1000 randomly selected)
x_II = Kd5_II ./ Kd6_II .* wp_II;
y_II = Kd7_II ./ Kd8_II .* wc_II;

x_UI = Kd5_UI ./ Kd6_UI .* wp_UI;
y_UI = Kd7_UI ./ Kd8_UI .* wc_UI;

figure(3001)
scatter(log10(x_II(1:1000)), log10(y_II(1:1000)), 10, 'filled')
legendInfo{1} = ['MI'];
hold on
scatter(log10(x_UI(1:1000)), log10(y_UI(1:1000)), 10, 'filled', 'MarkerFaceAlpha',0.8)
legendInfo{2} = ['J/U'];
xlim([-4,3])
xticks([-4:1:3])
ylim([-4,3])
yticks([-4:1:3])
xlabel("Log10 (Kd5/Kd6*Wp)")
ylabel("Log10 (Kd7/Kd8*Wc)")
pbaspect([1 1 1])
set(gca,'FontSize',18)
box off
x = [-4:3];
y = [-4:3];
plot(x,y, 'black', 'HandleVisibility', 'off')
legend(legendInfo,'location','northwest')
legend boxon  
title("MI vs. J/U")

% Fig. 9B: Histograms for parameter conditions differentiating monotonic increasing vs. U-shape curves
composite_II = (Kd7_II ./ Kd8_II .* wc_II) ./ (Kd5_II ./ Kd6_II .* wp_II);
composite_UI = (Kd7_UI ./ Kd8_UI .* wc_UI) ./ (Kd5_UI ./ Kd6_UI .* wp_UI);
figure (3002)
hold on
histogram(log10(composite_II), 'BinWidth', 0.5, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(composite_UI), 'BinWidth', 0.5, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(composite_II), 'BinWidth', 0.5, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
histogram(log10(composite_UI), 'BinWidth', 0.5, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
xlim([-5,5])
xticks([-4:1:4])
xlabel("Log10 ((Kd7/Kd8*Wc) / (Kd5/Kd6*Wp))")
ylabel('Frequency')
pbaspect([1 1 1])
set(gca,'FontSize',18)
legend(legendInfo,'location','northwest')
legend boxon  
title("Increasing vs. J-shape")


% Fig. 9C: U-shape curve magnitude vs Log10 ((Kd7/Kd8*Wc) / (Kd5/Kd6*Wp))
MagL = EE_curve_U(:,2) - min(EE_curve_U(:,2:end),[], 2);
MagL_normalized = (EE_curve_U(:,2) - min(EE_curve_U(:,2:end),[], 2)) ./ EE_curve_U(:,2);

MagR = EE_curve_U(:,end) - min(EE_curve_U(:,2:end),[], 2);
MagR_normalized = (EE_curve_U(:,end) - min(EE_curve_U(:,2:end),[], 2)) ./ EE_curve_U(:,end);

figure(3003)
scatter(log10(composite_UI), MagL_normalized, 10, 'filled')
xlabel("Log10 ((Kd7/Kd8*Wc) / (Kd5/Kd6*Wp))")
ylabel('Magnitude of U')
xlim([-2,4])
xticks([-2:1:4])
pbaspect([1 1 1])
set(gca,'FontSize',18)


% Fig. 9D: Scatter plot for parameter conditions differentiating monotonic decreasing vs. Bell-shape curves
figure(3004)
%scatter(log10(Kd7_DI ./ Kd5_DI .* wc_DI), log10(Kd6_DI ./ Kd8_DI .* wp_DI), 10, 'filled')
scatter(log10(Kd8_DI ./ Kd5_DI .* wc_DI), log10(Kd6_DI ./ Kd7_DI .* wp_DI), 10, 'filled')
legendInfo{1} = ['MD'];
hold on
%scatter(log10(Kd7_BI ./ Kd5_BI .* wc_BI), log10(Kd6_BI ./ Kd8_BI .* wp_BI), 10, 'filled', 'MarkerFaceAlpha',0.8)
scatter(log10(Kd8_BI ./ Kd5_BI .* wc_BI), log10(Kd6_BI ./ Kd7_BI .* wp_BI), 10, 'filled', 'MarkerFaceAlpha',0.8)
legendInfo{2} = ['Bell'];
xlim([-5,3])
xticks([-5:1:3])
ylim([-5,3])
yticks([-5:1:3])
xlabel("Log10 (Kd8/Kd5*Wc)")
ylabel("Log10 (Kd6/Kd7*Wp)")
pbaspect([1 1 1])
set(gca,'FontSize',18)
box off
x = [-5:3];
y = [-5:3];
plot(x,y, 'black', 'HandleVisibility', 'off')
legend(legendInfo,'location','northwest')
legend boxon  
title("MD vs. Bell")


% Fig. 9E: Histograms for parameter conditions differentiating monotonic decreasing vs. Bell-shape curves
composite_DI = (Kd6_DI ./ Kd7_DI .* wp_DI) ./ (Kd8_DI ./ Kd5_DI .* wc_DI);
composite_BI = (Kd6_BI ./ Kd7_BI .* wp_BI) ./ (Kd8_BI ./ Kd5_BI .* wc_BI);
figure (3005)
hold on
histogram(log10(composite_DI), 'BinWidth', 0.5, 'EdgeColor', 'none', 'DisplayName','MD')
histogram(log10(composite_BI), 'BinWidth', 0.5, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
histogram(log10(composite_DI), 'BinWidth', 0.5, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
histogram(log10(composite_BI), 'BinWidth', 0.5, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
xlim([-6,5])
xticks([-6:2:4])
ylim([0,130])
yticks([0:30:120])
xlabel("Log10 ((Kd6/Kd7*Wp) / (Kd8/Kd5*Wc))")
ylabel('Frequency')
pbaspect([1 1 1])
set(gca,'FontSize',18)
legend(legendInfo,'location','northwest')
legend boxon  
title("Decreasing vs. Bell-shape")


% Fig. 9F: Bell-shape curve magnitude vs Log10 ((Kd6/Kd7*Wp) / (Kd8/Kd5*Wc))
MagL = max(EE_curve_Bell(:,2:end),[], 2) - EE_curve_Bell(:,2);
MagL_normalized = (max(EE_curve_Bell(:,2:end),[], 2) - EE_curve_Bell(:,2)) ./ EE_curve_Bell(:,2);

MagR = max(EE_curve_Bell(:,2:end),[], 2) - EE_curve_Bell(:,end);
MagR_normalized = (max(EE_curve_Bell(:,2:end),[], 2) - EE_curve_Bell(:,end)) ./ EE_curve_Bell(:,end);

figure(3006)
scatter(log10(composite_BI), log10(MagL_normalized), 10, 'filled')
xlabel("Log10 ((Kd6/Kd7*Wp) / (Kd8/Kd5*Wc))")
ylabel('Magnitude of Bell)')
xlim([-2,5])
xticks([-2:1:5])
pbaspect([1 1 1])
set(gca,'FontSize',18)




% ------------Histograms for Wc and Wp for differentiating MI vs MD curves---------------%
% Fig. S10E: Histograms for Wp for differentiating MI vs MD curves
figure (30043)
hold on
histogram((wp_II), 'BinWidth', 0.02, 'EdgeColor', 'none', 'DisplayName','MI')
histogram((wp_DI), 'BinWidth', 0.02, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram((wp_II), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
histogram((wp_DI), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
xticks([0, 0.1, 0.5, 1])
xlabel("Wp")
ylabel('Frequency')
pbaspect([1 1 1])
set(gca,'FontSize',18)
pbaspect([1.4 1 1])
legend('location','northwest')
legend boxon  
title("MI vs. MD")

% Fig. S10F: Histograms for Wc for differentiating MI vs MD curves
figure (30042)
hold on
histogram((wc_II), 'BinWidth', 0.02, 'EdgeColor', 'none', 'DisplayName','MI')
histogram((wc_DI), 'BinWidth', 0.02, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram((wc_II), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
histogram((wc_DI), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
xticks([0, 0.1, 0.5, 1])
xlabel("Wc")
ylabel('Frequency')
pbaspect([1 1 1])
set(gca,'FontSize',18)
pbaspect([1.4 1 1])
legend('location','northwest')
legend boxon  
title("MI vs. MD")



% ------------Histograms for Wc and Wp for differentiating U-shape vs Bell-shape curves---------------%
% Fig. S10G: Histograms for Wp for differentiating U-shape vs Bell-shape curves
figure (30040)
hold on
histogram((wp_UI), 'BinWidth', 0.02, 'EdgeColor', 'none', 'DisplayName','J/U')
histogram((wp_BI), 'BinWidth', 0.02, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
histogram((wp_UI), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
histogram((wp_BI), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
xticks([0, 0.1, 0.5, 1])
xlabel("Wp")
ylabel('Frequency')
pbaspect([1 1 1])
set(gca,'FontSize',18)
pbaspect([1.4 1 1])
legend('location','northwest')
legend boxon  
title("J/U-shape vs. Bell-shape")

% Fig. S10H: Histograms for Wc for differentiating U-shape vs Bell-shape curves
figure (30041)
hold on
histogram((wc_UI), 'BinWidth', 0.02, 'EdgeColor', 'none', 'DisplayName','J/U')
histogram((wc_BI), 'BinWidth', 0.02, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
histogram((wc_UI), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'HandleVisibility', 'off')
histogram((wc_BI), 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'HandleVisibility', 'off')
xticks([0, 0.1, 0.5, 1])
xlabel("Wc")
ylabel('Frequency')
pbaspect([1 1 1])
set(gca,'FontSize',18)
pbaspect([1.4 1 1])
legend('location','northwest')
legend boxon  
title("J/U-shape vs. Bell-shape")





% ----------Histogram of individual parameters for the 4 curve shapes------------------%

% Fig. 10A
default_Kd5 = default_param.kb5/default_param.kf5;
% Kd5 - MI and J/U
figure(3005)
hold on
histogram(log10(Kd5_II/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(Kd5_UI/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(Kd5_II/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd5_UI/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (Kd5 fold change)")
ylabel('Frequency (x100)')
ylim([0,500])
yticks([0:250:500])
yticklabels({'0', '2.5', '5'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off


% Kd5 - MD and Bell
figure(30051)
hold on
histogram(log10(Kd5_DI/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd5_BI/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd5_DI/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(Kd5_BI/default_Kd5), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (Kd5 fold change)")
ylabel('Frequency (x100)')
yticks([0:50:100])
yticklabels({'0', '0.5', '1'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off


% Fig. 10B
default_Kd6 = default_param.kb6/default_param.kf6;
% Kd6 - MI and J/U
figure(3006)
hold on
histogram(log10(Kd6_II/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(Kd6_UI/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(Kd6_II/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd6_UI/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (Kd6 fold change)")
ylabel('Frequency (x100)')
ylim([0,500])
yticks([0:250:500])
yticklabels({'0', '2.5', '5'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off


% Kd6 - MD and Bell
figure(30061)
hold on
histogram(log10(Kd6_DI/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd6_BI/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd6_DI/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(Kd6_BI/default_Kd6), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (Kd6 fold change)")
ylabel('Frequency (x100)')
ylim([0,100])
yticks([0:50:100])
yticklabels({'0', '0.05', '0.1'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff  
legend off


% Fig. 10C
default_Kd7 = default_param.kb7/default_param.kf7;
% Kd7 - MI and J/U
figure(3007)
hold on
histogram(log10(Kd7_II/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(Kd7_UI/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(Kd7_II/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd7_UI/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (Kd7 fold change)")
ylabel('Frequency (x100)')
ylim([0,500])
yticks([0:250:500])
yticklabels({'0', '2.5', '5'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off

% Kd7 - MD and Bell
figure(30071)
hold on
histogram(log10(Kd7_DI/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd7_BI/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd7_DI/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(Kd7_BI/default_Kd7), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (Kd7 fold change)")
ylabel('Frequency (x100)')
ylim([0,100])
yticks([0:50:100])
yticklabels({'0', '0.5', '1'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff  
legend off


% Fig. 10D
default_Kd8 = default_param.kb8/default_param.kf8;
% Kd8 - MI and J/U
figure(3008)
hold on
histogram(log10(Kd8_II/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(Kd8_UI/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(Kd8_II/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd8_UI/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (Kd8 fold change)")
ylabel('Frequency (x100)')
ylim([0,500])
yticks([0:250:500])
yticklabels({'0', '2.5', '5'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off

% Kd8 - MD and Bell
figure(30081)
hold on
histogram(log10(Kd8_DI/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd8_BI/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd8_DI/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(Kd8_BI/default_Kd8), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (Kd8 fold change)")
ylabel('Frequency (x100)')
ylim([0,100])
yticks([0:50:100])
yticklabels({'0', '0.5', '1'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff  
legend off


% Fig. 10E
default_wp = default_param.wp;
% wp - MI and J/U
figure(3009)
hold on
histogram(wp_II/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(wp_UI/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(wp_II/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(wp_UI/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("wp")
ylabel('Frequency (x100)')
ylim([0,500])
yticks([0:250:500])
yticklabels({'0', '2.5', '5'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off

% wp - MD and Bell
figure(30091)
hold on
histogram(wp_DI/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(wp_BI/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(wp_DI/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(wp_BI/default_wp, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("wp")
ylabel('Frequency (x100)')
ylim([0,300])
yticks([0:100:300])
yticklabels({'0', '1', '2','3'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff  
legend off


% Fig. 10F
default_wc = default_param.wc;
% wc - MI and J/U
figure(3010)
hold on
histogram(wc_II/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(wc_UI/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(wc_II/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(wc_UI/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("wc")
ylabel('Frequency (x100)')
ylim([0,500])
yticks([0:250:500])
yticklabels({'0', '2.5', '5'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off


% wc - MD and Bell
figure(30101)
hold on
histogram(wc_DI/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(wc_BI/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(wc_DI/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(wc_BI/default_wc, 'Normalization','count', 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("wc")
ylabel('Frequency (x100)')
ylim([0,200])
yticks([0:100:200])
yticklabels({'0', '1', '2'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff  
legend off


% Fig. 10G
% k1 - MI and J/U
figure(3011)
hold on
histogram(log10(k1_II/default_param.k1), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(k1_UI/default_param.k1), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(k1_II/default_param.k1), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(k1_UI/default_param.k1), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (k1 fold change)")
ylabel('Frequency (x100)')
ylim([0,700])
yticks([0:350:700])
yticklabels({'0', '3.5', '7'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off

% k1 - MD and Bell
figure(30111)
hold on
histogram(log10(k1_DI/default_param.k1), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(k1_BI/default_param.k1), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(k1_DI/default_param.k1), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(k1_BI/default_param.k1), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (k1 fold change)")
ylabel('Frequency (x100)')
ylim([0,110])
yticks([0:50:100])
yticklabels({'0', '0.5', '1'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff  
legend off


% Fig. 10H
% k30 - MI and J/U
figure(3012)
hold on
histogram(log10(k30_II/default_param.k30), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(k30_UI/default_param.k30), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(k30_II/default_param.k30), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(k30_UI/default_param.k30), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (k30 fold change)")
ylabel('Frequency (x100)')
ylim([0,600])
yticks([0:300:600])
yticklabels({'0', '3', '6'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off

% k30 - MD and Bell
figure(30121)
hold on
histogram(log10(k30_DI/default_param.k30), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(k30_BI/default_param.k30), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(k30_DI/default_param.k30), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(k30_BI/default_param.k30), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (k30 fold change)")
ylabel('Frequency (x100)')
ylim([0,120])
yticks([0:50:100])
yticklabels({'0', '0.5', '1'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff  
legend off


% Fig. 10I
% k3 - MI and J/U
figure(3013)
hold on
histogram(log10(k3_II/default_param.k3), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(k3_UI/default_param.k3), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(k3_II/default_param.k3), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(k3_UI/default_param.k3), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (k3 fold change)")
ylabel('Frequency (x100)')
ylim([0,500])
yticks([0:250:500])
yticklabels({'0', '2.5', '5'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off

% k3 - MD and Bell
figure(30131)
hold on
histogram(log10(k3_DI/default_param.k3), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(k3_BI/default_param.k3), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(k3_DI/default_param.k3), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(k3_BI/default_param.k3), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (k3 fold change)")
ylabel('Frequency (x100)')
ylim([0,100])
yticks([0:50:100])
yticklabels({'0', '0.5', '1'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff  
legend off


% Fig. 10J
% Kd3 - MI and J/U
figure(3015)
hold on
histogram(log10(Kd3_II/default_param.Kd3), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(Kd3_UI/default_param.Kd3), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(Kd3_II/default_param.Kd3), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd3_UI/default_param.Kd3), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (Kd3 fold change)")
ylabel('Frequency (x100)')
ylim([0,300])
yticks([0:150:300])
yticklabels({'0', '1.5', '3'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off

% Kd3 - MD and Bell
figure(30151)
hold on
histogram(log10(Kd3_DI/default_param.Kd3), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd3_BI/default_param.Kd3), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(Kd3_DI/default_param.Kd3), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(Kd3_BI/default_param.Kd3), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (Kd3 fold change)")
ylabel('Frequency (x100)')
ylim([0,50])
yticks([0:25:50])
yticklabels({'0', '0.25', '0.5'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff  
legend off


% Fig. 10K
% n3 - MI and J/U
figure(3016)
hold on
histogram((n3_II/default_param.n3), 'Normalization','count', 'BinWidth', 0.02, 'EdgeColor', 'none', 'DisplayName','MI')
histogram((n3_UI/default_param.n3), 'Normalization','count', 'BinWidth', 0.02, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram((n3_II/default_param.n3), 'Normalization','count', 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram((n3_UI/default_param.n3), 'Normalization','count', 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("n3 fold change")
ylabel('Frequency (x100)')
ylim([0,400])
yticks([0:200:400])
yticklabels({'0', '2', '4'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off

% n3 - MD and Bell
figure(30161)
hold on
histogram((n3_DI/default_param.n3), 'Normalization','count', 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram((n3_BI/default_param.n3), 'Normalization','count', 'BinWidth', 0.02, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram((n3_DI/default_param.n3), 'Normalization','count', 'BinWidth', 0.02, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram((n3_BI/default_param.n3), 'Normalization','count', 'BinWidth', 0.02, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("n3 fold change")
ylabel('Frequency (x100)')
ylim([0,60])
yticks([0:25:50])
yticklabels({'0', '0.25', '0.5'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff  
legend off


% Fig. 10L
% CRtot - MI and J/U
figure(3017)
hold on
histogram(log10(CRtot_II/default_param.CRtot), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'DisplayName','MI')
histogram(log10(CRtot_UI/default_param.CRtot), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','J/U')
histogram(log10(CRtot_II/default_param.CRtot), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(CRtot_UI/default_param.CRtot), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
xlabel("Log10 (CRtot fold change)")
ylabel('Frequency (x100)')
ylim([0,300])
yticks([0:150:300])
yticklabels({'0', '1.5', '3'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff
legend off

% CRtot - MD and Bell
figure(30171)
hold on
histogram(log10(CRtot_DI/default_param.CRtot), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black',  'LineWidth', 0.75, 'DisplayName','')
histogram(log10(CRtot_BI/default_param.CRtot), 'Normalization','count', 'BinWidth', 0.1, 'DisplayStyle','stairs', 'EdgeColor', 'black', 'LineWidth', 0.75, 'DisplayName','')
histogram(log10(CRtot_DI/default_param.CRtot), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','MD')
histogram(log10(CRtot_BI/default_param.CRtot), 'Normalization','count', 'BinWidth', 0.1, 'EdgeColor', 'none', 'FaceAlpha',0.6, 'DisplayName','Bell')
xlabel("Log10 (CRtot fold change)")
ylabel('Frequency (x100)')
ylim([0,50])
yticks([0:25:50])
yticklabels({'0', '0.25', '0.5'})
pbaspect([1 0.7 1])
set(gca,'FontSize',32)
legend boxoff  
legend off


