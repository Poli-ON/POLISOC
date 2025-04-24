function dydt = ECMRC1order(time,y,t,I,TAU,soc,Q)
%Job of the function: Define the following differential equation:
%
% Variables:
% y(1) = i_R1
% y(2) = soc
%
% Equazione:
% dy(1)/dt=-1/R1C1+1/R1C1*i(t)
% dy(2)/dt=-i(t)/3600*Q

dydt = zeros(2,1);
I_int = interp1(t,I,time); % Interpolate the current I at time "time"

if length(TAU)>1 % In case you have defined a SOC dependent TAU
    TAU_int = interp1(soc,TAU,y(2));
else
    TAU_int=TAU;
end

dydt(1) =  -TAU_int.*y(1) + TAU_int.*I_int;
dydt(2) = -I_int./(3600*Q);

end