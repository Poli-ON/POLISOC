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
 
function Data = init_nh(v0,d0,T0,SigmaX0,SigmaV,SigmaW,param,ObType)

  % Initial state description: 
  % Save in the variable SOC0 the initial SOC computed by inverting the OCV
  % (assuming the cell is in equilibrium condition). The same job can be
  % obtained inverting the deformation in NMC and LCO cells (monotone THK
  % function).

  cyc=param.cyc;

  switch ObType
      case 1
          % Voltage inversion
          % SOC0  = interp1(param.OCV,param.SOC,min(max(v0,min(param.OCV)),max(param.OCV)));Data.zkInd = 1;
          % Deformation inversion
          SOC0  = interp1((param.DthkC{cyc}+flip(param.DthkD{cyc}))./2,param.SOC,min(max(d0,min((param.DthkC{cyc}+flip(param.DthkD{cyc}))./2)),max((param.DthkC{cyc}+flip(param.DthkD{cyc}))./2))); Data.zkInd = 1; % If d0 is lower/greater than the min/max bound of the registered SOC-THK relationship, takes the boundary value. (Avoid NaN with iterp1).
          Data.xhat  = SOC0; % initial state
      case 2
          ir0   = 0;                                                                        Data.irInd = 1;
          hk0   = 0;                                                                        Data.hekInd = 2;
          SOC0  = interp1(param.OCV,param.SOC,min(max(v0,min(param.OCV)),max(param.OCV)));  Data.zkInd = 3; % If v0 is lower/greater than the min/max bound of the registered SOC-OCV relationship, takes the boundary value. (Avoid NaN with iterp1).
          Data.xhat  = [ir0 hk0 SOC0]'; % initial state   
      case 3
          ir0   = 0;                                                                        Data.irInd = 1;
          hk0   = 0;                                                                        Data.hekInd = 2;
          % Voltage inversion
          % SOC0  = interp1(param.OCV,param.SOC,min(max(v0,min(param.OCV)),max(param.OCV)));Data.zkInd = 3;
          % Voltage + Deformation inversion
          SOC0  = (interp1(param.OCV,param.SOC,min(max(v0,min(param.OCV)),max(param.OCV))) + interp1((param.DthkC{cyc}+flip(param.DthkD{cyc}))./2,param.SOC,min(max(d0,min((param.DthkC{cyc}+flip(param.DthkD{cyc}))./2)),max((param.DthkC{cyc}+flip(param.DthkD{cyc}))./2))))/2;   Data.zkInd = 3; % If v0 or d0 is lower/greater than the min/max bound of the registered SOC-OCV or SOC-THK relationship, takes the boundary value. (Avoid NaN with iterp1). 
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
