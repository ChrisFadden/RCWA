function [Wref, Sref] = calcReflectionSide(Kx,Ky,KzR,N,W0,V0,grid)
  
  x1 = 1:N;
  x2 = N+1:2*N;
  y1 = x1;
  y2 = x2;

  Q = zeros(2*N,2*N);
  
  Q(x1,y1) = Kx * Ky;
  Q(x1,y2) = grid.erR*grid.urR*eye(N) - Kx.^2; %Check this...
  Q(x2,y1) = Ky.^2 - grid.erR*grid.urR*eye(N);
  Q(x2,y2) = -Ky * Kx;
  
  Q = Q / grid.urR;

  Wref = zeros(2*N,2*N);
  Wref(x1,y1) = eye(N);
  Wref(x2,y2) = eye(N);
  
  %The -1 are fudge factors
  %kind of makes sense for going opposite direction...
  lam = W0;
  lam(x1,y1) = -1*(sqrt(-1)*KzR);
  lam(x2,y2) = -1*(sqrt(-1)*KzR);
      
  Vref = Q * inv(lam);
  
  A = inv(W0)*Wref + inv(V0)*Vref;
  B = inv(W0)*Wref - inv(V0)*Vref;

  Sref.S11 = -inv(A)*B;
  Sref.S12 = 2*inv(A);
  Sref.S21 = 0.5*(A - B*inv(A)*B);
  Sref.S22 = B*inv(A); 
end
