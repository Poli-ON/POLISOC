% Job of the function:
%    run EFK filter for SOC estimation with multiple observers
%
% Inputs:
%   ObType: Type of observer to be used (1: Deformation, 2: Voltage, 3:
%   Hybrid)
%   param: Set of cell parameters
%   Meas: Structure containing measurement data of the test to estimate SOC
%   BatType: Type of battery (0: LCO, 1: LFP, 2: NMC, 3: NMC Michigan dataset)
%
% Output:
%   sochat: SOC estimation
%   socbound: 3-sigma estimation bounds
%   OutModel: Model Output 


function [sochat,socbound,OutModel]=run_nh(ObType,param,Meas,BatType)

%% Covariance values
% Set the Covariance values
switch BatType
    case 0 %LCO
        switch ObType
            case 1
                SigmaX0 = diag(0.0025); % uncertainty of initial state (SOC)
                SigmaV  = 1e-3; % Output equation: Uncertainty of deformation sensor
                SigmaW  = 1e-2; % State equation: Uncertainty of current sensor
            case 2
                SigmaX0 = diag([0.001 1e-8 0.0025]); % uncertainty of initial state (current resistor, hysteresis, SOC)
                SigmaV  = 1e-1; % Output equation: Uncertainty of voltage sensor
                SigmaW  = 1e-2; % State equation: uncertainty of current sensor
            case 3
                SigmaX0 = diag([0.001 1e-8 0.0025]); % uncertainty of initial state
                SigmaV  = diag([1e-3; 1e-1]); %  Output equation: Uncertainty of deformation and voltage sensor
                SigmaW  = 1e-2; % State equation: Uncertainty of current sensor
        end
    case 1 %LFP
        switch ObType
            case 1
                SigmaX0 = diag(0.0025); % uncertainty of initial state (SOC)
                SigmaV  = 1e-3; % Output equation: Uncertainty of deformation sensor
                SigmaW  = 1e-2; % State equation: Uncertainty of current sensor
            case 2
                SigmaX0 = diag([0.001 1e-8 0.0025]); % uncertainty of initial state (current resistor, hysteresis, SOC)
                SigmaV  = 0.5e-2; % Output equation: Uncertainty of voltage sensor
                SigmaW  = 1e-2; % State equation: uncertainty of current sensor
            case 3
                SigmaX0 = diag([0.001 1e-8 0.0025]); % uncertainty of initial state
                SigmaV  = diag([1e-3; 0.5e-2]); %  Output equation: Uncertainty of deformation and voltage sensor
                SigmaW  = 1e-2; % State equation: Uncertainty of current sensor
        end
    case 2 %NMC
        switch ObType
            case 1
                SigmaX0 = diag(0.0025); % uncertainty of initial state (SOC)
                SigmaV  = 1e-3; % Output equation: Uncertainty of deformation sensor
                SigmaW  = 1e-2; % State equation: Uncertainty of current sensor
            case 2
                SigmaX0 = diag([0.001 1e-8 0.0025]); % uncertainty of initial state (current resistor, hysteresis, SOC)
                SigmaV  = 1e-3; % Output equation: Uncertainty of voltage sensor
                SigmaW  = 1e-2; % State equation: uncertainty of current sensor
            case 3
                SigmaX0 = diag([0.001 1e-8 0.0025]); % uncertainty of initial state
                SigmaV  = diag([1e-3; 1e-2]); %  Output equation: Uncertainty of deformation and voltage sensor
                SigmaW  = 1e-2; % State equation: Uncertainty of current sensor
        end
    case 3 %NMC Michigan dataset
        switch ObType
            case 1
                SigmaX0 = diag(0.0025); % uncertainty of initial state (SOC)
                SigmaV  = 1e-3; % Output equation: Uncertainty of deformation sensor
                SigmaW  = 1e-2; % State equation: Uncertainty of current sensor
            case 2
                SigmaX0 = diag([0.001 1e-8 0.0025]); % uncertainty of initial state (current resistor, hysteresis, SOC)
                SigmaV  = 1e-3; % Output equation: Uncertainty of voltage sensor
                SigmaW  = 1e-2; % State equation: uncertainty of current sensor
            case 3
                SigmaX0 = diag([0.001 1e-8 0.0025]); % uncertainty of initial state
                SigmaV  = diag([1e-3; 1e-3]); %  Output equation: Uncertainty of deformation and voltage sensor
                SigmaW  = 1e-2; % State equation: Uncertainty of current sensor
        end
end

% Time variables
n = length(Meas.Time);
deltat = Meas.Time(3)-Meas.Time(2);           %Sampling interval [s]

% Allocate variables for results
sochat = zeros(n,1);
if ObType == 3
    OutModel = zeros(n,2);
else
    OutModel = zeros(n,1);
end
socbound = zeros(n,1);

% Fixed temperature
T=20;

% Create the Data structure and initialize variables using first
% voltage or deformation measurement
Data = init_nh(Meas.Voltage(1),Meas.Deformation(1),T,SigmaX0,SigmaV,SigmaW,param,ObType);

% Time loop of the filter, updated each time step "k" with a new measurement
% value.
hwait = waitbar(0,'Computing...'); 
for k = 1:length(Meas.Time)
  dk = Meas.Deformation(k);   % "measure" deformation
  vk = Meas.Voltage(k);       % "measure" voltage
  ik = Meas.Current(k);       % "measure" current
  Tk = T;                     % "measure" temperature
  
 % Compute state variables according to the algorith mode:
  switch ObType
      case 1
          [sochat(k),socbound(k),OutModel(k),Data] = iter_deformation_nh(dk,ik,Tk,deltat,Data);
      case 2
          [sochat(k),socbound(k),OutModel(k),Data] = iter_voltage_nh(vk,ik,Tk,deltat,Data);
      case 3
          [sochat(k),socbound(k),OutModel(k,:),Data] = iter_hybrid_nh(dk,vk,ik,Tk,deltat,Data);
  end

 % update waitbar periodically (slow procedure)
  if mod(k,1000)==0
    waitbar(k/length(Meas.Current),hwait);
  end
end
close(hwait);



end
