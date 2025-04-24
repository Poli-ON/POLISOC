clear variables


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
        load param_LFP11.mat
        param=param_LFP11;
        param.cyc=1;
    case 2 %NMC
        t0=14;                   % Nominal battery thickness
        load param_NMC1.mat
        param=param_NMC1;
        param.cyc=1;
end

cyc = param.cyc;
Q  = param.Q;
Gm  = param.Gm;
Mm  = param.Mm;

%% Load test data

switch BatType
    case 0 %LCO
        load("Meas_LCO_DST1.mat")
    case 1 % LFP
        % load("Meas_LFP11_DST1.mat")
        load("Meas_LFP11_DriveCycle_1.mat")
    case 2 % NMC
        load("Meas_NMC1_DST1.mat")
        % load("Meas_NMC2_DriveCycle_2.mat")
end

I=Meas.TrueCurrent;

% Retrieve deformation at SOC 100%
maxthk=param.DthkC{param.cyc}(end);
% Tolgo contributo termico
Meas.Deformation=Meas.Deformation-param.alfa*t0.*(Meas.Temperature-20);
% Faccio partire il modello dalla deformazione corrispondente a SOC 100%
Meas.Deformation = Meas.Deformation - Meas.Deformation(1) + maxthk;


%% MODEL

dt=diff(Meas.Time);
SOC=ones(length(Meas.Time),1);
for i=2:length(Meas.Time)
    SOC(i)=SOC(i-1)-I(i).*dt(i-1)./(3600*Q);
end

y0=[1;0];
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[time_t,y] = ode45(@(time,y) EMModefun(time,y,Meas.Time,I,Gm,Q),Meas.Time,y0,opts);

% Deformation without hysteresis
tnh=DTHKfromSOC(y(:,1),param,param.cyc);
% Evaluate the magnitude of the deformation hysteresis envelope at the present SOC
Mmk=interp1(param.SOC,param.Mm,y(:,1)); 
indmin=y(:,1)<0.01;
indmax=y(:,1)>0.99;
Mmk(indmin)=param.Mm(1);
Mmk(indmax)=param.Mm(end);
% Deformation with hysteresis
t=tnh+y(:,2).*Mmk;


rmse=sqrt(sum((Meas.Deformation-t).^2)./length(t))

createfigure(Meas.Time, [Meas.Deformation,t],time_t, y(:,1),rmse)
% createfigure(Meas.SOC, [Meas.Deformation,t],SOC, y(:,1))

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
ylabel('Deformation [mm]','FontName','Times New Roman');

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
str = ['Rmse = ', sprintf('%.3f', rmse), ' mm'];
annotation(figure1,'textbox',...
    [0.158928571428571 0.179952379499164 0.420000009587833 0.0847619062151228],...
    'String',str,'FitBoxToText','on','FontSize',17,...
    'FontName','Times New Roman',...
    'BackgroundColor',[1 1 1]);
end



