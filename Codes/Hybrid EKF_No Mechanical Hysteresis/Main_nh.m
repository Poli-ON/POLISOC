clear variables

%% Load battery type parameter database

BatType = 0;% 0 LCO, 1 LFP, 2 NMC
switch BatType
    case 0 %LCO
        t0=9.2;                  % Nominal battery thickness
        load param_LCO.mat
        param=param_LCO;
        param.cyc=44;            % Set the aging level 1: first RPT, 44: last RPT. Reference: 10.5281/zenodo.14914172
    case 1 %LFP
        t0=27;                   % Nominal battery thickness
        load param_LFP11.mat
        param=param_LFP11;
        param.cyc=1;
    case 2 %NMC
        t0=14;                   % Nominal battery thickness
        load param_NMC2.mat
        param=param_NMC2;
        param.cyc=1;
end


%% Load test data

switch BatType
    case 0 %LCO
        load("Meas_LCO_DST2.mat")
    case 1 % LFP
        % load("Meas_LFP11_DST2.mat")
        load("Meas_LFP11_DriveCycle_2.mat")
    case 2 % NMC
        load("Meas_NMC1_DST1.mat")
        % load("Meas_NMC2_DriveCycle_1.mat")
end

% Add a random error and a bias to the current
Meas.Current = Meas.TrueCurrent + Meas.TrueCurrent.*(rand(length(Meas.Current),1)-0.5)./10 + sign(Meas.Current).*max(Meas.Current)./50;

% Compensate temperature on deformation measurements
Meas.Deformation = Meas.Deformation - t0.*param.alfa*(Meas.Temperature-Meas.Temperature(1));
Meas.Deformation=Meas.Deformation-Meas.Deformation(end);

n = length(Meas.Time);
deltat = Meas.Time(3)-Meas.Time(2);           %Sampling interval [s]


%% Compute the reference "TRUE" SOC with Coulomb counting for comparison

Ah = zeros(n,1);
Ah_current_sensor = zeros(n,1);
for i=2:n
    Ah(i) = Ah(i-1)+Meas.TrueCurrent(i).*(deltat/3600);
    Ah_current_sensor(i) = Ah_current_sensor(i-1)+Meas.Current(i).*(deltat/3600);
end


switch BatType
    case 0 %LCO
        nominalCap                      = Ah(end); 
    case 1 %LFP
        nominalCap                      = Ah(end); 
    case 2 %NMC
        nominalCap                      = Ah(end); 
end
SOC0_cc                         = 1;
soc                             = (SOC0_cc * nominalCap - Ah)./nominalCap;                  % True soc computed with coulomb counting
soc_current                     = (SOC0_cc * nominalCap - Ah_current_sensor)./nominalCap;   % Fictious SOC considering fictitious sensor errors
T                               = 20;                                                       % Fixed temperarure

%% RUN

% The first field of the function "run" controls how to run the algorithm:
% 1: Run the algorith in Deformation mode
% 2: Run the algorith in Voltage mode
% 3: Run the algorith in Hybrid mode

% In this case, run the algorithm in all the three modes
for i=1:3
    [sochat{i},socbound{i},OutModel{i}]=run_nh(i,param,Meas,BatType);
    err{i}=sqrt(mean((100*(soc-sochat{i})).^2));
end

%% PLOT
color=[0.00,0.45,0.74; 0.85,0.33,0.10; 0.93,0.69,0.13];

% Plot SOC estimation
figure();
p1=subplot(1,1,1);
plot(Meas.Time/60,100*soc,'k','linewidth',1); hold on;
for i=1:3
    plot(Meas.Time/60,100*sochat{i},'linewidth',1,'Color',color(i,:))
    hold on
end
% Possibly plot also the fictious coulomb counting SOC, affected by current
% sensor synthetic errors.
% plot(Meas.Time/60,100*soc_current,'--k'); hold on;
set(p1,'TickLabelInterpreter', 'tex','FontSize',17,'FontName','Times New Roman')
xlabel('Time [min]'); ylabel('SOC [%]');
legend('Reference',['Deformation. Error= ' num2str(err{1},'%0.2f') '%'],['Voltage. Error= ' num2str(err{2},'%0.2f') '%'],['Hybrid. Error= ' num2str(err{3},'%0.2f') '%']); ylim([0 120]); grid on
ylim([0 100])



% Plot  model
figure();
p1=subplot(2,1,1);
plot(Meas.Time/60,OutModel{1},'linewidth',0.75); hold on; plot(Meas.Time/60,OutModel{3}(:,1),'linewidth',0.75); hold on; plot(Meas.Time/60,Meas.Deformation,'k','linewidth',0.75);
set(p1,'TickLabelInterpreter', 'tex','FontSize',17,'FontName','Times New Roman')
xlabel('Time [min]'); 
ylabel('Deformation [mm]')
grid on

p2=subplot(2,1,2);
plot(Meas.Time/60,OutModel{2},'linewidth',0.75,'Color',color(1,:)); hold on; plot(Meas.Time/60,OutModel{3}(:,2),'linewidth',0.75,'Color',color(2,:)); hold on; plot(Meas.Time/60,Meas.Voltage,'k','linewidth',0.75);
legend('Model - Single output','Model - Double output','Measurements','Location','south')
set(p2,'TickLabelInterpreter', 'tex','FontSize',17,'FontName','Times New Roman')
xlabel('Time [min]');
ylabel('Voltage [V]')
% ylim([2.5 4.2])
ylim([2.4 3.5])
grid on


% Plot estimation error and bounds
figure(); 
p1=subplot(1,1,1);
for i=1:3
    plot(Meas.Time/60,100*(soc-sochat{i}),'Color',color(i,:),'linewidth',1); hold on
    h = plot([Meas.Time/60; NaN; Meas.Time/60],[100*socbound{i}; NaN; -100*socbound{i}],'--','Color',color(i,:),'linewidth',1);
end
xlabel('Time [min]'); ylabel('SOC error [%]'); ylim([-10 10]); 
set(p1,'TickLabelInterpreter', 'tex','FontSize',17,'FontName','Times New Roman')
legend('Error','Bounds','location','northwest'); 
grid on


% Print info in the command window

% % Display RMS estimation error to command window
% fprintf('RMS SOC estimation error = %g%%\n',err);

% Display bounds errors to command window
% ind = find(abs(soc-sochat)>socbound);
% fprintf('Percent of time error outside bounds = %g%%\n',...
%         length(ind)/length(soc)*100);

