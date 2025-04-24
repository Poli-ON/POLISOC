function docv=dOCVfromSOC(soc,param)

% This function computes an approximation to the derivative of the
% ocv relationship (docv) of a cell with respect to state of charge.
%
% Inputs:  soc = scalar or matrix of cell state of charge points
%        Fun = data structure produced by processDynamic
% Output: docv = scalar or vector of derivatives of open circuit voltage

  soccol = soc(:); % force soc to be col-vector
  SOC = param.SOC(:); % force to be col vector... 
  OCV0 = param.OCV(:); % force to be col vector... 
  dOCV0 = gradient(OCV0)./gradient(SOC);

  diffSOC=SOC(2)-SOC(1); % spacing between OCV points - assume uniform
  docv=zeros(size(soccol)); % initialize output to zero
  indL=find(soccol <= SOC(1)); % indices of socs below model-stored data
  indH=find(soccol >= SOC(end)); % and of socs above model-stored data
  indR=find(soccol > SOC(1) & soccol < SOC(end)); % the rest of them
  indN=isnan(soccol); % if input is "not a number" for any locations

  % For V < lowest SoC point, extrapolate using first table values 
  if ~isempty(indL)
    dv = (dOCV0(2)) - (dOCV0(1));
    docv(indL)= (soccol(indL)-SOC(1)).*dv/diffSOC + dOCV0(1); 
  end

  % For V > highest SoC point, extrapolate using last table values
  if ~isempty(indH)
    dv = (dOCV0(end)) - (dOCV0(end-1));
    docv(indH) = (soccol(indH)-SOC(end)).*dv/diffSOC + dOCV0(end);
  end

  % for normal soc range, manually interpolate (faster than "interp1")
  I4=(soccol(indR)-SOC(1))/diffSOC; % using linear interpolation
  I5=floor(I4); I45 = I4-I5; omI45 = 1-I45;
  docv(indR)=dOCV0(I5+1).*omI45 + dOCV0(I5+2).*I45;
  docv(indN)=0; % replace NaN SOCs with zero voltage
  docv = reshape(docv,size(soc)); % output is same shape as input