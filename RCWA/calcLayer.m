function Si = calcLayer(Kx,Ky,N,W0,V0,grid,device,layer)
  
  %Get parameters for current layer
  ERC = device.ERC(:,:,layer);
  URC = device.URC(:,:,layer);
    
  %Setup matrix indexing
  x1 = 1:N;
  x2 = N+1:2*N;
  y1 = x1;
  y2 = x2;

  Q = zeros(2*N,2*N);
  P = zeros(2*N,2*N);
  
  %Compute PQ Matrices
  P(x1,y1) = Kx * inv(ERC) * Ky;
  P(x1,y2) = URC - Kx * inv(ERC) * Kx;
  P(x2,y1) = Ky * inv(ERC)*Ky - URC;
  P(x2,y2) = - Ky * inv(ERC) * Kx;
  
  Q(x1,y1) = Kx * inv(URC) * Ky;
  Q(x1,y2) = ERC - Kx * inv(URC) * Kx;
  Q(x2,y1) = Ky * inv(URC) * Ky - ERC;
  Q(x2,y2) = - Ky * inv(URC) * Kx;
    
  %Compute Eigenmodes of layer
  %Eigenvalues /vectors can be out of order, but as long as they are
  %together, the end Scattering Matrices will be the same.
  OM2 = P*Q;
  [W, LAM2] = eig(OM2);

  LAM = sqrt(LAM2);
  V = Q*W*inv(LAM);

  
  %Compute S-matrix temp variables
  A = inv(W) * W0 + inv(V)*V0;
  B = inv(W) * W0 - inv(V)*V0;
  X = expm(-LAM * (2*pi / grid.lam0) * grid.L(layer));
  
  Si.S11 = inv(A - X*B*inv(A)*X*B)*(X*B*inv(A)*X*A - B);
  Si.S12 = inv(A - X*B*inv(A)*X*B)*X*(A - B*inv(A)*B);
  Si.S21 = Si.S12;
  Si.S22 = Si.S11;  
end
