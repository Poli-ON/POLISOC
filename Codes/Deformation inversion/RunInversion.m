clear variables

%% Load battery type parameter database

BatType = 2;% 0 LCO, 1 LFP, 2 NMC, 3 NMC Michigan dataset

switch BatType
    case 0 %LCO
        t0=9.2;             % Nominal battery thickness
        load param_LCO.mat
        param=param_LCO;
        param.cyc=1;       % Set the aging level 1: first RPT, 44: last RPT. Reference: 10.5281/zenodo.14914172
    case 1 %LFP
        t0=27;              % Nominal battery thickness
        load param_LFP10.mat
        param=param_LFP10;
        param.cyc=1;
    case 2 %NMC
        t0=14;              % Nominal battery thickness
        load param_NMC1.mat
        param=param_NMC1;
        % param.DthkD{1}=param.DthkD{1}./1.05; % Attivare per prove a 40°C
        % param.DthkC{1}=param.DthkC{1}./1.05; % Attivare per prove a 40°C
        param.cyc=1;
    case 3 % NMC Michigan dataset
        t0=4;                   % Nominal battery thickness
        load param_NMC19.mat
        param=param_NMC19;
        param.cyc=6;            % 6 aging steps: 1: fresh, 2: 106 cycles, 3: 213 c. 4: 319 c. 5: 424 c.  6: 449c. (EOL)
end


%% Load test data

switch BatType
    case 0 %LCO
        load("Meas_LCO_DST2.mat")
    case 1 % LFP
        load("Meas_LFP10_DST1.mat")
        % load("Meas_LFP10_DriveCycle_1.mat")
    case 2 % NMC
        % load("Meas_NMC1_DST1.mat")
        load("Meas_NMC1_DriveCycle_1.mat")
    case 3 % NMC Michigan dataset
        load("Meas_NMC19.mat")
        Meas=Meas{6};
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
    case 3 %NMC Michigan dataset
        nominalCap                      = 2*Ah(end); 
end
SOC0_cc0                        = 1;
soc_cc                          = (SOC0_cc0 * nominalCap - Ah)./nominalCap;                  % True soc computed with coulomb counting
soc_current                     = (SOC0_cc0 * nominalCap - Ah_current_sensor)./nominalCap;   % Fictious SOC considering fictitious sensor errors
T                               = 20;


%% RUN
tic
Meas.Deformation(Meas.Deformation>max(param.DthkD{param.cyc}))=max(param.DthkD{param.cyc}); % Force deformation greater than the SOC-deformation characteristics to the max of the characteristic (SOC = 1)
Meas.Deformation(Meas.Deformation<min(param.DthkD{param.cyc}))=min(param.DthkD{param.cyc}); % Force deformation lower than the SOC-deformation characteristics to the min of the characteristic (SOC = 0)
soc=interp1((flip(param.DthkD{param.cyc}))/1,param.SOC,Meas.Deformation); % Invert the deformation characteristics recorded during discharge
% soc= soc - soc(1) + SOC0_cc0; % Force SOC0_cc0 to be initial soc.
toc

err=sqrt(mean((100*(soc-soc_cc)).^2));

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
