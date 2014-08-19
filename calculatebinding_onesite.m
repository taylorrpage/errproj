function [d, a, da] = calculatebinding_onesite(dAdded, aAdded, Ka)
%[d, a, da] = calculatebinding_onesite(dAdded, aAdded, Ka)
%Returns the equilibrium values, d, a, and da for a one site binding model given
%the analytical concentrations of D, dAdded, and A, aAdded, and a binding
%constant, Ka
%
%Concentrations, dAdded and aAdded, should be in units of uM.
%Ka should be in units of 1/M

Kd = 1/Ka * 10^6;

sum = dAdded + aAdded + Kd;
da  = (sum - sqrt(sum^2 - 4*dAdded*aAdded))/2;

d   = dAdded - da;
a   = aAdded - da;