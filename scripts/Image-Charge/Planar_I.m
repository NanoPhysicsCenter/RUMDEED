function [I] = Planar_I(V_0, w_theta, R, d)
%PLANAR_I Summary of this function goes here
%   Detailed explanation goes here

F = V_0/d;
FN_st = FN_J(w_theta, F);

A = R^2*pi;
I = A*FN_st;
end

