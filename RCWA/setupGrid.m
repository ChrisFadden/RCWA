function [g] = setupGrid( fileName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fb = importdata(fileName);
g.lam0 = fb(1);
g.Lx = fb(2);
g.Ly = fb(3);
g.Nx = fb(4);
g.Ny = fb(5);
g.dx = g.Lx / g.Nx;
g.dy = g.Ly / g.Ny;
g.urR = fb(6);
g.erR = fb(7);
g.urT = fb(8);
g.erT = fb(9);
L = fb(10);

Ldim = [];
for i = 1:L
    Ldim = [fb(10 + i) Ldim];
end
g.L = Ldim;

end

