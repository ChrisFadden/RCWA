function [W0, V0, SG] = calcFreeSpace(Kx,Ky,N)
  
  Kz = conj(sqrt(eye(N) - Kx.^2 - Ky.^2));
    
  x1 = 1:N;
  x2 = N+1:2*N;
  y1 = x1;
  y2 = x2;

  Q = zeros(2*N,2*N);
  
  Q(x1,y1) = Kx * Ky;
  Q(x1,y2) = eye(N) - Kx.^2; 
  Q(x2,y1) = Ky.^2 - eye(N);
  Q(x2,y2) = -Kx * Ky;
    
  W0 = zeros(2*N,2*N);
  W0(x1,y1) = eye(N);
  W0(x2,y2) = eye(N);
  
  LAM = W0;
  LAM(x1,y1) = (sqrt(-1)*Kz);
  LAM(x2,y2) = (sqrt(-1)*Kz);
    
  V0 = Q*inv(LAM);
  
  SG.S11 = zeros(2*N,2*N);
  SG.S21 = eye(2*N);
  SG.S12 = eye(2*N);
  SG.S22 = zeros(2*N);
end
