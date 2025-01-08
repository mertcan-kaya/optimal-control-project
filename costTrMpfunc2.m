function y = costTrMpfunc2(G)

A = [0,1;-1,0];
B = [0;1];

Q = [G(1),0;0,G(2)];
P = G(3);

G
R = icare(A,B,Q,P,[],[],[])

K = R*B/P;

k1 = K(1);
k2 = K(2);

omegan = sqrt((k1 + 1));
zeta = k2/(2*omegan);

omegad = omegan*sqrt(1-zeta^2);
phi = acos(zeta);

tr = real((pi-phi)/omegad);

tp = real(pi/omegad);
Mp = exp(-zeta*omegan*tp)*(-cosh(sqrt(zeta^2-1)*omegan*tp)+sinh(sqrt(zeta^2-1)*omegan*tp)*(zeta-omegan)/(sqrt(zeta^2-1)));

y(1,1) = abs(tr)
y(2,1) = abs(Mp)