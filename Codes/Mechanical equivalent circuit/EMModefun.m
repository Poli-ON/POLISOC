function dydt = EMModefun(time,y,t,I,G,Q)

% Variables:
% y(1) = soc
% y(2) = hysteresis
%
% Equation:
% dy(1)/dt=-i(t)/3600*Q
% dy(2)/dt=-abs(I_int.*G./Q)*h + abs(I_int.*G./Q)*sgn(i(t))

dydt = zeros(2,1);
I_int = interp1(t,I,time); % Interpolate the current I at time "time"

ab  = abs(I_int.*G./Q);

dydt(1) = -I_int./(3600*Q);
dydt(2) =  -ab.*y(2) + ab.*sign(I_int);

end