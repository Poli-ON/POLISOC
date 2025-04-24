function dydt = ECMRC2order(time,y,t,I,TAU1,TAU2,Q)
%Job of the function: Define the following differential equation:
%
% Variables:
% y(1) = i_R1
% y(2) = iR2
% y(3) = soc
%
% Equazione:
% dy(1)/dt=-1/R1C1+1/R1C1*i(t)
% dy(2)/dt=-1/R2C2+1/R2C2*i(t)
% dy(3)/dt=-i(t)/3600*Q

dydt = zeros(3,1);
I_int = interp1(t,I,time); % Interpolate the current I at time "time"
dydt(1) =  -TAU1.*y(1) + TAU1.*I_int;
dydt(2) =  -TAU2.*y(2) + TAU2.*I_int;
dydt(3) = -I_int./(3600*Q);

end
