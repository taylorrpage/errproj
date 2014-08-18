function Y = expfun(param, t)

    y0 = param(1);
    A0 = param(2);
    k1 = param(3);
    c  = param(4);

    %Check if first time point is early enough for k1.
    [t, hasAddedT] = checkTi(k1, t);
    
    %Numerical integration
    [~, yOut] = ode45(@diffFun, t, y0);
    Y = A0 .* yOut + c;
    
    %Remove first point if preliminary ti used.
    if tcheck
        Y = Y(2:end);
    end
    
    %System of diff eq
    function dy = diffFun(~, y)
        dy = -k1*y(1);
    end


end