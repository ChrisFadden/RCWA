function [Wtrn, Strn] = calcTransmissionSide(Kx,Ky,KzT,N,W0,V0,grid)
  
  x1 = 1:N;
  x2 = N+1:2*N;
  y1 = x1;
  y2 = x2;

  Q = zeros(2*N,2*N);
  
  Q(x1,y1) = Kx * Ky;
  Q(x1,y2) = grid.erT*grid.urT*eye(N) - Kx.^2; %Check this...
  Q(x2,y1) = Ky.^2 - grid.erT*grid.urT*eye(N);
  Q(x2,y2) = -Ky * Kx;
  
  Q = Q / grid.urT;

  Wtrn = zeros(2*N,2*N);
  Wtrn(x1,y1) = eye(N);
  Wtrn(x2,y2) = eye(N);

  lam = W0;
  lam(x1,y1) = sqrt(-1)*KzT;
  lam(x2,y2) = sqrt(-1)*KzT;
    
  Vtrn = Q*inv(lam);
  
  A = inv(W0)*Wtrn + inv(V0)*Vtrn;
  B = inv(W0)*Wtrn - inv(V0)*Vtrn;

  Strn.S11 = B*inv(A);
  Strn.S12 = 0.5 * (A - B*inv(A)*B);
  Strn.S21 = 2*inv(A);
  Strn.S22 = -inv(A)*B; 
end
