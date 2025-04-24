% Job of the function:
%    Initialize a data structure, used by the extended Kalman
%    filter to store its own state and associated data.
%
% Inputs:
%   v0: Initial cell voltage
%   d0: Initial cell deformation (thickness)
%   T0: Initial cell temperature
%   SigmaX0: Initial state uncertainty covariance matrix
%   SigmaV: Covariance of measurement noise
%   SigmaW: Covariance of process noise
%
% Output:
%   Data: Data structure used by EKF algorithm code
 
function Data = init(v0,d0,T0,SigmaX0,SigmaV,SigmaW,param,ObType)

  % Initial state description: 
  % Save in the variable SOC0 the initial SOC computed by inverting the OCV
  % (assuming the cell is in equilibrium condition). The same job can be
  % obtained inverting the deformation in NMC and LCO cells (monotone THK
  % function).

  switch ObType
      case 1
          % Voltage inversion
          SOC0  = interp1(param.OCV,param.SOC,v0);                                  Data.zkInd = 1;
          % Deformation inversion
          % SOC0  = interp1((param.DthkC{1}+flip(param.DthkD{1}))./2,param.SOC,d0); Data.zkInd = 1;
          Data.xhat  = SOC0; % initial state
      case 2
          ir0   = 0;                                                                Data.irInd = 1;
          hk0   = 0;                                                                Data.hkInd = 2;
          SOC0  = interp1(param.OCV,param.SOC,v0);                                  Data.zkInd = 3;
          Data.xhat  = [ir0 hk0 SOC0]'; % initial state   
      case 3
          ir0   = 0;                                                                Data.irInd = 1;
          hk0   = 0;                                                                Data.hkInd = 2;
          % Voltage inversion
          SOC0  = interp1(param.OCV,param.SOC,v0);                                  Data.zkInd = 3;
          % Deformation inversion
          % SOC0  = interp1((param.DthkC{1}+flip(param.DthkD{1}))./2,param.SOC,d0);   Data.zkInd = 3;
          Data.xhat  = [ir0 hk0 SOC0]'; % initial state
  end

  % Covariance values
  Data.SigmaX = SigmaX0;
  Data.SigmaV = SigmaV;
  Data.SigmaW = SigmaW;
  Data.Qbump = 10;
  
  % previous value of current
  Data.priorI = 0;
  Data.signIk = 0;
  
  % store model data structure too
  Data.param = param;

end