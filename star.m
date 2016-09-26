function S = star(SA,SB)

sizeM = size(SA.S12);
sizeM = sizeM(1);

D = SA.S12 * inv(eye(sizeM) - SB.S11*SA.S22);
F = SB.S21 * inv(eye(sizeM) - SA.S22*SB.S11);

S.S11 = SA.S11 + D*SB.S11*SA.S21;
S.S12 = D*SB.S12;
S.S21 = F*SA.S21;
S.S22 = SB.S22 + F*SA.S22*SB.S12;

end
