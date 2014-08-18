function [t, tcheck] = ticheck(k1, t)

tcheck = false;
ti = log(1/k1) - 2;

if t(1) > 10^ti
    t = [10^ti; t];
    tcheck = true;
end