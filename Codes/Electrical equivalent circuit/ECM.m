clear variables

ECMorder = 1;

%% Load battery type parameter database

BatType = 2;% 0 LCO, 1 LFP, 2 NMC
switch BatType
    case 0 %LCO
        t0=9.2;                  % Nominal battery thickness
        load param_LCO.mat
        param=param_LCO;
        param.cyc=44;            % Set the aging level. 1: first RPT, 44: last RPT. Reference: 10.5281/zenodo.14914172
    case 1 %LFP
        t0=27;                   % Nominal battery thickness
        load param_LFP10.mat
        param=param_LFP10;
        param.cyc=1;
    case 2 %NMC
        t0=14;                   % Nominal battery thickness
        load param_NMC2.mat
        param=param_NMC2;
        param.cyc=1;
end

cyc = param.cyc;
Q  = param.Q;
Ge  = param.Ge;
Me  = param.Me;
M0 = param.M0;
R1 = param.R1;
C1 = param.C1;
R0 = param.R0;
eta = 1;
TAU1 = 1/(R1*C1);
soc=[];


%% Load test data

switch BatType
    case 0 %LCO
        load("Meas_LCO_DST1.mat")
    case 1 % LFP
        % load("Meas_LFP10_DST2.mat")
        load("Meas_LFP10_DriveCycle_1.mat")
    case 2 % NMC
        % load("Meas_NMC2_DST1.mat")
        load("Meas_NMC2_DriveCycle_2.mat")
end

n = length(Meas.Time);
deltat = Meas.Time(3)-Meas.Time(2);           %Sampling interval [s]


%% MODEL
switch ECMorder
    case 1
        y0=[0;interp1(param.OCV,param.SOC,Meas.Voltage(1))];
        opts = odeset('RelTol',1e-4,'AbsTol',1e-6);
        [time_t,y] = ode45(@(time,y) ECMRC1order(time,y,Meas.Time,Meas.TrueCurrent,TAU1,soc,Q),Meas.Time,y0,opts);
        v=OCVfromSOC(y(:,2),param)-R1.*y(:,1)-R0*Meas.TrueCurrent;
    case 2
        y0=interp1(param.OCV,param.SOC,Meas.Voltage(1));
        opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [time_t,y] = ode45(@(time,y) ECMRC2order(time,y,Meas.Time,Meas.TrueCurrent,TAU1,TAU2,Q),Meas.Time,y0,opts);
        v=OCVfromSOC(y(:,2),param)-R1.*y(:,1)-R0*Meas.TrueCurrent;
end


rmse=sqrt(sum((Meas.Voltage-v).^2)./length(v))

createfigure(Meas.Time, [Meas.Voltage,v], time_t, y(:,2),rmse)

function createfigure(X1, YMatrix1, X2, Y1,rmse)
%CREATEFIGURE(X1, YMatrix1, X2, Y1)
%  X1:  vector of plot x data
%  YMATRIX1:  matrix of plot y data
%  X2:  vector of plot x data
%  Y1:  vector of plot y data

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
colororder([0.929 0.694 0.125]);

% Activate the left side of the axes
yyaxis(axes1,'left');
% Create multiple line objects using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',1);
set(plot1(1),'DisplayName','Measurement','Color',[0 0 0],'LineStyle','-');
set(plot1(2),'DisplayName','Model','Color',[0.85,0.33,0.10],'LineStyle','-');

% Create ylabel
ylabel('Voltage [V]','FontName','Times New Roman');
% ylim([2.75 4.2])
ylim([2.5 3.65])

% Set the remaining axes properties
set(axes1,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(axes1,'right');
% Create plot
plot(X2,Y1,'DisplayName','SOC','LineWidth',1.5,...
    'Color',[0 0.447058823529412 0.741176470588235],'LineStyle',':');
% Create ylabel
ylabel('State of Charge [-]','FontName','Times New Roman',...
    'Color',[0 0.447058823529412 0.741176470588235]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0 1]);

% Set the remaining axes properties
set(axes1,'YColor',[0 0.447058823529412 0.741176470588235]);
% Create xlabel
xlabel('Time [s]','FontName','Times New Roman');

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontName','Times New Roman','FontSize',17);
% Create legend
legend(axes1,'show');

% Create textbox
str = ['Model rmse = ', sprintf('%.3f', rmse), ' V'];
annotation(figure1,'textbox',...
    [0.158928571428571 0.179952379499164 0.420000009587833 0.0847619062151228],...
    'String',str,'FitBoxToText','on','FontSize',17,...
    'FontName','Times New Roman',...
    'BackgroundColor',[1 1 1]);
end





