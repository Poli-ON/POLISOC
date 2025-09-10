% Job of the function:
%    Perform one iteration of the deformation-based extended Kalman filter using the new
%    measured data.
%
% Inputs:
%   dk: Present measured (noisy) cell deformation
%   ik: Present measured (noisy) cell current
%   Tk: Present temperature
%   deltat: Sampling interval
%   Data: Data structure initialized by init, updated by iter
%
% Output:
%   zk: SOC estimate for this time sample
%   zkbnd: 3-sigma estimation bounds
%   yat: Output deformation (model) 
%   Data: Data structure used to store persistent variables
 
function [zk,zkbnd,yhat,Data] = iter_deformation(dk,ik,Tk,deltat,Data)

  % Load the cell model parameters
  param=Data.param;
  cyc = param.cyc;

  Q  = param.Q;
  Gm  = param.Gm;
  Gm=0;
  Mm  = param.Mm;
  eta = 1;
  if ik<0, ik=ik*eta; end % Multiply by Coulomb efficiency when discharge (ik>0)
  
  % Get data stored in Data structure
  I = Data.priorI; %i_k+omega_k
  SigmaX = Data.SigmaX;
  SigmaV = Data.SigmaV;
  SigmaW = Data.SigmaW;
  xhat = Data.xhat;
  zkInd = Data.zkInd;
  hmkInd = Data.hmkInd;
  if abs(ik)>Q/100, Data.signIk = sign(ik); end % Update the current value if it exceeds C/100
  signIk = Data.signIk;
  
  % EKF Step 0: Compute Ahat[k-1], Bhat[k-1]
  nx = length(xhat); Ahat = zeros(nx,nx); Bhat = zeros(nx,1);
  Ahat(zkInd,zkInd) = 1; Bhat(zkInd) = -deltat/(3600*Q);
  Ah  = exp(-abs(I*Gm*deltat/(3600*Q)));  % hysteresis factor
  Ahat(hmkInd,hmkInd) = Ah;
  % Write B-matrix of the state equation
  B = [Bhat, 0*Bhat]; % write B matrix of the states
  Bhat(hmkInd) = -abs(Gm*deltat/(3600*Q))*Ah*(1+sign(I)*xhat(hmkInd));
  B(hmkInd,2) = Ah-1;
  
  % Step 1: State estimate time update
  xhat = Ahat*xhat + B*[I; sign(I)];
  % help mantein robustness, Ensure that no value falls outside the range [-0.05, 1.05]
  xhat(zkInd) = min(1.05,max(-0.05,xhat(zkInd)));

  % Step 2: Error covariance time update
  SigmaX = Ahat*SigmaX*Ahat' + Bhat*SigmaW*Bhat';
  
  % Step 3: Output estimate
  if ~any(Mm) % To speed up don't use interp1 if mechanical hysteresis is not considered (Mm= [0 ... 0]).
      Mmk=0;
  else
      Mmk=interp1(param.SOC,Mm,min(max(xhat(zkInd),0),1)); % Calculate amplitude of deformation hysteresis envelope at present soc, and assign 
  % the boundary envelope magnitude in case soc is outside the range [0 1]
  end
  yhat = DTHKfromSOC(xhat(zkInd),param,cyc) + Mmk*xhat(hmkInd);
  
  % Step 4: Compute the estimator gain matrix
  Chat = zeros(1,nx);
  Chat(zkInd) = dDTHKfromSOC(xhat(zkInd),param,cyc);
  Chat(hmkInd) = Mmk;
  Dhat = 1;
  SigmaY = Chat*SigmaX*Chat' + Dhat*SigmaV*Dhat';
  L = SigmaX*Chat'/SigmaY;
  
  % Step 5: State estimate measurement update
  r = dk - yhat;    % error between measured deformation and model
  if r^2 > 100*SigmaY, L(:)=0.0; end 
  xhat = xhat + L*r;
  xhat(hmkInd) = min(1,max(-1,xhat(hmkInd))); % Help maintain robustness, Ensure that no value falls outside the range [-1, 1]
  xhat(zkInd) = min(1.05,max(-0.05,xhat(zkInd)));% Help maintain robustness, Ensure that no value falls outside the range [-0.05, 1.05]
  
  % Step 6: Error covariance measurement update
  SigmaX = SigmaX - L*SigmaY*L';
  % Q-bump code 
  % If sigmaY is small, I trust the deformation measurement, so I want a large L.
  % If the innovation (error r) is large, there is high uncertainty in the state estimate (sigmaX).
  if r^2 > 4*SigmaY % bad deformation estimate by 2 std. devs, bump Q 
    fprintf('Bumping SigmaX\n');
    SigmaX(zkInd,zkInd) = SigmaX(zkInd,zkInd)*Data.Qbump;
  end
  [~,S,V] = svd(SigmaX);
  HH = V*S*V';
  SigmaX = (SigmaX + SigmaX' + HH + HH')/4; % Help maintain robustness
  
  % Save data in Data structure for next timestep
  Data.priorI = ik;
  Data.SigmaX = SigmaX;
  Data.xhat = xhat;
  zk = xhat(zkInd);
  zkbnd = 3*sqrt(SigmaX(zkInd,zkInd));
end
