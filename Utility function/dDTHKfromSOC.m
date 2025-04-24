function ddthk=dDTHKfromSOC(soc,param,cyc)

% This function computes an approximation to the derivative of the
% reversible deformation relationship (ddthk) of a cell with respect to 
% state of charge.
%
% Inputs:  soc = scalar or matrix of cell state of charge points
%        Fun = data structure 
% Output: dDthk = scalar or vector of derivatives of reversible

  soccol = soc(:); % force soc to be col-vector
  SOC = param.SOC(:); % force to be col vector... 

  % Compute the reference SOC-THK relationship as the mean curve between
  % the reversible deformation measured during a complete charge (DthkC) 
  % and discharge (DthkD)
  DTHK0 = (param.DthkC{cyc}(:)+flip(param.DthkD{cyc}(:)))./2; % force to be column vector

  dDTHK0 = gradient(DTHK0)./gradient(SOC);

  diffSOC=SOC(2)-SOC(1); % spacing between OCV points - assume uniform
  ddthk=zeros(size(soccol)); % initialize output to zero
  indL=find(soccol <= SOC(1)); % indices of socs below model-stored data
  indH=find(soccol >= SOC(end)); % and of socs above model-stored data
  indR=find(soccol > SOC(1) & soccol < SOC(end)); % the rest of them
  indN=isnan(soccol); % if input is "not a number" for any locations

  % For ddthk < lowest SoC point, extrapolate using first table values 
  if ~isempty(indL)
    dv = (dDTHK0(2)) - (dDTHK0(1));
    ddthk(indL)= (soccol(indL)-SOC(1)).*dv/diffSOC + dDTHK0(1);
  end

  % For ddthk > highest SoC point, extrapolate using last table values
  if ~isempty(indH)
    dv = (dDTHK0(end)) - (dDTHK0(end-1));
    ddthk(indH) = (soccol(indH)-SOC(end)).*dv/diffSOC + dDTHK0(end);
  end

  % for normal soc range, manually interpolate (faster than "interp1")
  I4=(soccol(indR)-SOC(1))/diffSOC; % using linear interpolation
  I5=floor(I4); I45 = I4-I5; omI45 = 1-I45;
  ddthk(indR)=dDTHK0(I5+1).*omI45 + dDTHK0(I5+2).*I45;
  ddthk(indN)=0; % replace NaN SOCs with zero deformation
  ddthk = reshape(ddthk,size(soc)); % output is same shape as input