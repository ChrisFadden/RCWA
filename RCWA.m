%*************************************************
% Simulation Parameters
%*************************************************
close all;
clc;
clear all;

%****************
% RCWA Parameters
%****************
Nx = 512;                     %x points in real-space grid
NH = 3;               %number of spatial harmonics for x, y
N = NH^2;             %size of matrices

%*******************
%   Read Data files
%*******************
[grid,Ny]   = setupGrid('grid.mat',Nx);
src = setupSrc('src.mat',N);
[device,grid] = setupDevice('device.mat',grid,Nx,Ny,NH);

%*******************************
% Calculate Kx, Ky, Kz Matrices 
%*******************************
[Kx, Ky, KzT, KzR] = calcK(src,grid,NH);

%******************
% Calc Eigenmodes
%******************
[W0,V0,SG] = calcFreeSpace(Kx,Ky,N);
[Wref, Sref] = calcReflectionSide(Kx,Ky,KzR,N,W0,V0,grid); 
[Wtrn, Strn] = calcTransmissionSide(Kx,Ky,KzT,N,W0,V0,grid); 

%Initialize Global Scattering Matrix
SG.S11 = zeros(2*N,2*N);
SG.S12 = eye(2*N);
SG.S21 = eye(2*N);
SG.S22 = zeros(2*N,2*N);

%********************************
%   Calc Global Scattering Matrix
%********************************
 for layer = 1:length(grid.L)                                                                                                                        
    Si = calcLayer(Kx,Ky,N,W0,V0,grid,device,layer);
    SG = star(SG,Si);
 end
 SG = star(Sref,SG);
 SG = star(SG,Strn);

%**********************
%   Post-Processing
%**********************

%Calculate Source Mode Coefficients
csrc = inv(Wref) * src.esrc;

%**********************
%   REFLECTED FIELDS
%**********************
eref = Wref * SG.S11 * csrc;
rx = eref(1:N);
ry = eref(N+1:end);
rz = - inv(KzR) * (Kx*rx + Ky*ry);
R = abs(rx).^2 + abs(ry).^2 + abs(rz).^2;
kzInc = cos(src.theta)*sqrt(grid.erR * grid.urR);

R = real(-KzR / kzInc) * R;
REF = sum(R)
R = reshape(R,[NH,NH]);

%***********************
%   TRANSMITTED FIELDS
%***********************
etrn = Wtrn * SG.S21 * csrc;
tx = etrn(1:N);
ty = etrn(N+1:end);
tz = -inv(KzT) * (Kx*tx + Ky*ty);
T = abs(tx).^2 + abs(ty).^2 + abs(tz).^2;

T = real((grid.urR / grid.urT) * KzT / kzInc) * T;
TRN = sum(T)
T = reshape(T,[NH,NH]);





















