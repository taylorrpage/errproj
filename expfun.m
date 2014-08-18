function Y = expfun(param, t)

    y0 = param(1);
    A0 = param(2);
    k1 = param(3);
    c  = param(4);

    %Check if first time point is early enough for k1.
    [t, tcheck] = ticheck(k1, t);
    
    %Numerical integration
    [~, yout] = ode45(@difeq, t, y0);
    Y = A0 .* yout + c;
    
    %Remove first point if preliminary ti used.
    if tcheck
        Y = Y(2:end);
    end
    
    %System of diff eq
    function dy = difeq(~, y)
        dy = -k1*y(1);
    end


end