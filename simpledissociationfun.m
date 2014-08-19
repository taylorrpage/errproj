function Y = simpledissociationfun(concentrations, param, t)
%t is a vertical vector
%param = [a0 kD kQ Ka kOff]
%concentrations = [zntAdded fe3Added]
%param units:
%   a0    complicated
%   kD    1/s
%   kQ    1/s
%   KaExp unitless - KaExp = log(Ka), with Ka units of 1/M
%   kOn   1/(M*s)

zntAdded  = concentrations(1);
fe3Added  = concentrations(2);
a0        = param(1);
kD        = param(2);
kQ        = param(3);
KaExp     = param(4);
kOn       = param(5);

Ka   = 10^KaExp;
kOff = kOn/Ka;
kOn  = kOn*10^-6; %kOn must be in uM-1s-1.

[znt0, fe30, zntfe30] = calculatebinding_onesite(zntAdded, fe3Added, Ka);

%Check if first time point is early enough for kQ.
[t, hasAddedT] = checkTi(kQ+kD, t); %zntfe3 decays as the sum of kQ and kD.

[~, concentrationsY] = ode45(@diffFun, t, [znt0 zntfe30]);

Y = a0 .* (concentrationsY(:, 1) + concentrationsY(:, 2));

%If an extra time point was added, the computed simulation point should be
%removed.
if hasAddedT
    Y = Y(2:end);
end

function dy = diffFun(~, y)
    fe3 = fe30 - y(2);
    dy = [-kD*y(1) - kOn*fe3*y(1) + kOff*y(2); ...
          -(kD + kQ + kOff)*y(2) + kOn*fe3*y(1)];
end

end