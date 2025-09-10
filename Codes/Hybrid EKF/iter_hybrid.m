% Job of the function:
%    Perform one iteration of the hybrid extended Kalman filter using the new
%    measured data.
%
% Inputs:
%   dk: Present measured (noisy) cell deformation (thickness
%   vk: Present measured (noisy) cell voltage
%   ik: Present measured (noisy) cell current
%   Tk: Present temperature
%   deltat: Sampling interval
%   Data: Data structure initialized by initEKF, updated by iterEKF
%
% Output:
%   zk: SOC estimate for this time sample
%   zkbnd: 3-sigma estimation bounds
%   yat: Output deformation (model) 
%   Data: Data structure used to store persistent variables
 
function [zk,zkbnd,yhat,Data] = iter_hybrid(dk,vk,ik,Tk,deltat,Data)

  % Load the cell model parameters
  param=Data.param;
  cyc = param.cyc;
  Q  = param.Q;
  Ge  = param.Ge;
  Gm  = param.Gm;
  Mm  = param.Mm;
  Me  = param.Me;
  M0 = param.M0;
  if iscell(param.R1)
      R0 = param.R0{cyc};
      R1 = param.R1{cyc};
      C1 = param.C1{cyc};
  else
      R0 = param.R0;
      R1 = param.R1;
      C1 = param.C1;
  end
  % R2 = param.R2(Tk);
  % C2 = param.C2(Tk);
  RC1= exp(-deltat./abs(R1*C1))';
  % RC2= exp(-deltat./abs(R2*C2))';
  RC = [RC1];
  eta = 1;
  if ik<0, ik=ik*eta; end % Multiply by Coulomb efficiency when discharge (ik>0)
  
  % Get data stored in Data structure
  I = Data.priorI; %i_k+omega_k
  SigmaX = Data.SigmaX;
  SigmaV = Data.SigmaV;
  SigmaW = Data.SigmaW;
  xhat = Data.xhat;
  irInd = Data.irInd;
  hekInd = Data.hekInd;
  zkInd = Data.zkInd;
  hmkInd = Data.hmkInd;
  if abs(ik)>Q/100, Data.signIk = sign(ik); end % Update the current value if it exceeds C/100
  signIk = Data.signIk;
  
  % EKF Step 0: Compute Ahat[k-1], Bhat[k-1]
  nx = length(xhat); Ahat = zeros(nx,nx); Bhat = zeros(nx,1);
  Ahat(zkInd,zkInd) = 1; Bhat(zkInd) = -deltat/(3600*Q);
  Ahat(irInd,irInd) = diag(RC); Bhat(irInd) = 1-RC(:); % May have more RC branches, Save values R_jC_j in a vector
  Ahe  = exp(-abs(I*Ge*deltat/(3600*Q)));  % hysteresis Voltage factor 
  Ahm  = exp(-abs(I*Gm*deltat/(3600*Q)));  % hysteresis deformation factor 
  Ahat(hekInd,hekInd) = Ahe;
  Ahat(hmkInd,hmkInd) = Ahm;
  B = [Bhat, 0*Bhat]; % write B matrix of the states
  Bhat(hekInd) = -abs(Ge*deltat/(3600*Q))*Ahe*(1+sign(I)*xhat(hekInd));
  Bhat(hmkInd) = -abs(Gm*deltat/(3600*Q))*Ahm*(1+sign(I)*xhat(hmkInd));
  B(hekInd,2) = Ahe-1;
  B(hmkInd,2) = Ahm-1;
  
  % Step 1: State estimate time update
  xhat = Ahat*xhat + B*[I; sign(I)]; 
  % help mantein robustness, 
  xhat(hekInd) = min(1,max(-1,xhat(hekInd))); % Ensure that no value falls outside the range  [-1,1]
  xhat(zkInd) = min(1.05,max(-0.05,xhat(zkInd)));% Ensure that no value falls outside the range  [-0.05,1.05]

  % Step 2: Error covariance time update
  SigmaX = Ahat*SigmaX*Ahat' + Bhat*SigmaW*Bhat';
  
  % Step 3: Output estimate
  if ~any(Mm) % To speed up don't use interp1 if mechanical hysteresis is not considered (Mm= [0 ... 0]).
      Mmk=0;
  else
  Mmk=interp1(param.SOC,Mm,min(max(xhat(zkInd),0),1)); % Calculate amplitude of deformation hysteresis envelope at present soc, and assign 
  % the boundary envelope magnitude in case soc is outside the range [0 1]
  end
  if ~any(Me) % To speed up don't use interp1 if electrical hysteresis is not considered (Me= [0 ... 0]).
      Mek=0;
  else
      Mek=interp1(param.SOC,Me,min(max(xhat(zkInd),0),1)); % Calculate amplitude of voltage hysteresis envelope at present soc, and assign 
  % the boundary envelope magnitude in case soc is outside the range [0 1]
  end
  nout = 2;
  yhat = zeros(nout,1);
  yhat(1) = DTHKfromSOC(xhat(zkInd),param,cyc) + Mmk*xhat(hmkInd);
  yhat(2) = OCVfromSOC(xhat(zkInd),param) + M0*signIk + ...
         Mek*xhat(hekInd) - R1*xhat(irInd) - R0*ik;
  
  % Step 4: Compute the estimator gain matrix
  Chat = zeros(nout,nx);
  Chat(1,zkInd) = dDTHKfromSOC(xhat(zkInd),param,cyc);
  Chat(1,hmkInd) = Mmk;
  Chat(2,zkInd) = dOCVfromSOC(xhat(zkInd),param);
  Chat(2,hekInd) = Mek;
  Chat(2,irInd) = -R1;
  Dhat = diag([1,1]);
  SigmaY = Chat*SigmaX*Chat' + Dhat*SigmaV*Dhat';
  L = SigmaX*Chat'/SigmaY;
  
  % Step 5: State estimate measurement update
  r = zeros(nout,1);
  r(1) = dk - yhat(1); % residual.  % error between measured deformatioj and model
  r(2) = vk - yhat(2); % residual.  % error between measured voltage and model
  if r.^2 > 100*SigmaY, L(:)=0.0; end 
  xhat = xhat + L*r; 
  xhat(hekInd) = min(1,max(-1,xhat(hekInd))); % Help maintain robustness, Ensure that no value falls outside the range [-1, 1]
  xhat(hmkInd) = min(1,max(-1,xhat(hmkInd))); % Help maintain robustness, Ensure that no value falls outside the range [-1, 1]
  xhat(zkInd) = min(1.05,max(-0.05,xhat(zkInd)));
  
  % Step 6: Error covariance measurement update
  SigmaX = SigmaX - L*SigmaY*L';
  % Q-bump code
  % If sigmaY is small, I trust the deformation measurement, so I want a large L.
  % If the innovation (error r) is large, there is high uncertainty in the state estimate (sigmaX).
  if r(1)^2 > 4*SigmaY(1,1) ||  r(2)^2 > 4*SigmaY(2,2) % bad model estimate by 2 std. devs, bump Q 
    fprintf('Bumping SigmaX\n');
    SigmaX(zkInd,zkInd) = SigmaX(zkInd,zkInd)*Data.Qbump;
  end
  [~,S,V] = svd(SigmaX);
  HH = V*S*V';
  SigmaX = (SigmaX + SigmaX' + HH + HH')/4; % Help maintain robustness
  
  % Save data in Data structure for next time step
  Data.priorI = ik;
  Data.SigmaX = SigmaX;
  Data.xhat = xhat;
  zk = xhat(zkInd);
  zkbnd = 3*sqrt(SigmaX(zkInd,zkInd));
  yhat=yhat';
end
