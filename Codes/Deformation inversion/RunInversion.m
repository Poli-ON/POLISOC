clear variables

%% Load battery type parameter database

BatType = 2;% 0 LCO, 1 LFP, 2 NMC

switch BatType
    case 0 %LCO
        t0=9.2;             % Nominal battery thickness
        load param_LCO.mat
        param=param_LCO;
        param.cyc=44;       % Set the aging level 1: first RPT, 44: last RPT. Reference: 10.5281/zenodo.14914172
    case 1 %LFP
        t0=27;              % Nominal battery thickness
        load param_LFP10.mat
        param=param_LFP10;
        param.cyc=1;
    case 2 %NMC
        t0=14;              % Nominal battery thickness
        load param_NMC2.mat
        param=param_NMC2;
        param.cyc=1;
end


%% Load test data

switch BatType
    case 0 %LCO
        load("Meas_LCO_DST1.mat")
    case 1 % LFP
        load("Meas_LFP10_DST1.mat")
        % load("Meas_LFP10_DriveCycle_1.mat")
    case 2 % NMC
        % load("Meas_NMC2_DST1.mat")
        load("Meas_NMC2_DriveCycle_3.mat")
end

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
        nominalCap                      = Ah(end); % Or, param.Q;
    case 1 %LFP
        nominalCap                      = Ah(end); % Or, param.Q;
    case 2 %NMC
        nominalCap                      = Ah(end); % Or, param.Q;
end
SOC0_cc                         = 1;
soc_cc                          = (SOC0_cc * nominalCap - Ah)./nominalCap;  % Calculate the SOC using Coloumb Counting for comparison
soc_current                     = (SOC0_cc * nominalCap - Ah_current_sensor)./nominalCap;
T                               = 20;


%% RUN

soc=interp1((flip(param.DthkD{param.cyc}))/1,param.SOC,Meas.Deformation); % Invert the deformation characteristics recorded during discharge
err=sqrt(mean((100*(soc(145:end)-soc_cc(145:end))).^2));


%% PLOT

figure();
p1=subplot(1,1,1);
plot(Meas.Time./60,soc_cc.*100,'k','linewidth',1); hold on; plot(Meas.Time./60,soc.*100,'linewidth',1)
set(p1,'TickLabelInterpreter', 'tex','FontSize',17,'FontName','Times New Roman')
xlabel('Time [min]')
ylabel('SOC [%]')
legend('Reference','Estimation')
grid on
ylim([0 100])