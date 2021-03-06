function [t, hasAddedT] = checkTi(k1, t)
%[t, hasAddedT] = checkTi(k1, t)
%Checks if the first time point is early enough for an accurate numerical
%solution given a decay constant.

hasAddedT = false;
ti = log(1/k1) - 2;

if t(1) > 10^ti
    t = [10^ti; t];
    hasAddedT = true;
end

end