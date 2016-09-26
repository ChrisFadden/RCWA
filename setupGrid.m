function [ g, Ny ] = setupGrid( fileName, Nx )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%*****************
% Layer Parameters
%*****************
g.lam0 = 2 * 10^-2;   %Free space wavelength
g.Lx = 17.5 * 10^(-3);  %period in x-direction
g.Ly = 15 * 10^(-3);    %period in y-direction
d1 = 5 * 10^(-3);     %thickness of layer 1
d2 = 3 * 10^(-3);     %thickness of layer 2

Ny = round(Nx * g.Ly / g.Lx); %y points in real-space grid

%***********************
%   Discrete Parameters
%***********************
g.dx = g.Lx / Nx;       %grid resolution in x
g.dy = g.Ly / Ny;       %grid resolution in y

g.L = [d1; d2];


end

