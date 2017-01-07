%*************************************************
% Simulation Parameters
%*************************************************
close all;
clc;
clear all;

%****************
% RCWA Parameters
%****************
NH = 3;               %number of spatial harmonics for x, y
if(mod(NH,2) == 0)
    error('Number of Harmonics must be odd');
end
N = NH^2;             %size of matrices

%*******************
%   Read Data files
%*******************
[grid] = setupGrid('./device/Triangle/Grid.dat');
src = setupSrc('src.mat',N);
device.ERC(N,N,length(grid.L)) = 0;
device.URC(N,N,length(grid.L)) = 0;
for layer = 1:length(grid.L)
    fn1 = './device/Triangle/';
    [device] = setupDevice(fn1,grid,NH,device,layer);
end

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
%   
% %**********************
% %   Post-Processing
% %**********************
% 
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
Gamma = rx + ry + rz;
kzInc = cos(src.theta)*sqrt(grid.erR * grid.urR);
R = real(-KzR / kzInc) * R;
Gamma = reshape(Gamma,[NH,NH]);
REF = norm(sum(sum(R)));


%***********************
%   TRANSMITTED FIELDS
%***********************
etrn = Wtrn * SG.S21 * csrc;
tx = etrn(1:N);
ty = etrn(N+1:end);
tz = -inv(KzT) * (Kx*tx + Ky*ty);
Tau = tx + ty + tz; 
T = abs(tx).^2 + abs(ty).^2 + abs(tz).^2;
T = real((grid.urR / grid.urT) * KzT/kzInc) * T;
Tau = reshape(Tau,[NH,NH]);
TRN = norm(sum(sum(T)));

%Gamma + Tau
 
%Plotting
Sx = reshape(tx,[NH,NH]);
normFactor = sum(src.esrc);
Ex = zeros(grid.Nx,grid.Ny);
ErTest = zeros(grid.Nx,grid.Ny);
ERCtest = device.ERC(:,:,1);
Me = (-(NH-1)/2):(NH-1)/2;
Ne = (-(NH-1)/2):(NH-1)/2;
for x = 1:grid.Nx
    for y = 1:grid.Ny
        for m = 1:NH
            for n = 1:NH
                kx = 2*pi*Me(m) / grid.Lx;
                ky = 2*pi*Ne(n) / grid.Ly;
                xr = x*grid.dx;
                yr = y*grid.dy;
                Ex(x,y) = Ex(x,y) + Sx(m,n)*exp(-1j*(kx*xr + ky*yr));
            end
        end
    end
end

%imagesc(real(Ex) / normFactor) 
%colorbar;











