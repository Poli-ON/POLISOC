function dthk=DTHKfromSOC(soc,param,cyc)

% Computes the reversible thickness (dthk) of the battery for a given state 
% of charge
%
% Inputs: soc = scalar or matrix of cell state of charge points
%       param = data structure containing ocv and reversible deformation
% Output: ocv = scalar or vector of open circuit voltages


  soccol = soc(:); % force soc to be column vector
  SOC = param.SOC(:); % force to be column vector

  % Compute the reference SOC-THK relationship as the mean curve between
  % the reversible deformation measured during a complete charge (DthkC) 
  % and discharge (DthkD)
  DTHK0 = (param.DthkC{cyc}(:)+flip(param.DthkD{cyc}(:)))./2; % force to be column vector

  diffSOC=SOC(2)-SOC(1); % spacing between SOC points - assume uniform
  dthk=zeros(size(soccol)); % initialize output to zero
  indL=find(soccol <= SOC(1)); % indices of socs below model-stored data
  indH=find(soccol >= SOC(end)); % and of socs above model-stored data
  indR=find(soccol > SOC(1) & soccol < SOC(end)); % the rest of them
  indN=isnan(soccol); % if input is "not a number" for any locations

  % For dthk < lowest SoC point, extrapolate using first table values
  if ~isempty(indL)
    dv = DTHK0(2) - DTHK0(1);
    dthk(indL)= (soccol(indL)-SOC(1)).*dv/diffSOC + DTHK0(1);
  end

  % For dthk > highest SoC point, extrapolate using last table values
  if ~isempty(indH)
    dv = DTHK0(end) - DTHK0(end-1);
    dthk(indH) = (soccol(indH)-SOC(end)).*dv/diffSOC + DTHK0(end);
  end

  % for normal soc range, manually interpolate (faster than "interp1")
  I4=(soccol(indR)-SOC(1))/diffSOC; % using linear interpolation
  I5=floor(I4); I45 = I4-I5; omI45 = 1-I45;
  dthk(indR)=DTHK0(I5+1).*omI45 + DTHK0(I5+2).*I45;
  dthk(indN)=0; % replace NaN SOCs with zero d
  dthk = reshape(dthk,size(soc)); % output is same shape as input
