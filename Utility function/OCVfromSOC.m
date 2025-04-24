function ocv=OCVfromSOC(soc,param)

% Computes the open circuit voltage of the battery for a given state of 
% charge
%
% Inputs: soc = scalar or matrix of cell state of charge points
%       param = data structure containing ocv and reversible deformation
% Output: ocv = scalar or vector of open circuit voltages

  soccol = soc(:); % force soc to be column vector
  SOC = param.SOC(:); % force to be column vector
  OCV0 = param.OCV(:); % force to be column vector


  diffSOC=SOC(2)-SOC(1); % spacing between SOC points - assume uniform
  ocv=zeros(size(soccol)); % initialize output to zero
  indL=find(soccol <= SOC(1)); % indices of socs below model-stored data
  indH=find(soccol >= SOC(end)); % and of socs above model-stored data
  indR=find(soccol > SOC(1) & soccol < SOC(end)); % the rest of them
  indN=isnan(soccol); % if input is "not a number" for any locations

  % For V < lowest SoC point, extrapolate using first table values
  if ~isempty(indL)
    dv = OCV0(2) - OCV0(1);
    ocv(indL)= (soccol(indL)-SOC(1)).*dv/diffSOC + OCV0(1);
  end

  % For V > highest SoC point, extrapolate using last table values
  if ~isempty(indH)
    dv = OCV0(end) - OCV0(end-1);
    ocv(indH) = (soccol(indH)-SOC(end)).*dv/diffSOC + OCV0(end);
  end

  % for normal soc range, manually interpolate (faster than "interp1")
  I4=(soccol(indR)-SOC(1))/diffSOC; % using linear interpolation
  I5=floor(I4); I45 = I4-I5; omI45 = 1-I45;
  ocv(indR)=OCV0(I5+1).*omI45 + OCV0(I5+2).*I45;
  ocv(indN)=0; % replace NaN SOCs with zero voltage
  ocv = reshape(ocv,size(soc)); % make output same shape as input
