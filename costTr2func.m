function y = costTr2func(G)

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

tshft = -acosh((omegan-zeta)/sqrt(1-2*zeta*omegan+omegan^2))/(sqrt(zeta^2-1)*omegan);

tr = tshft+real((pi-phi)/omegad);

% tp = tshft+real(pi/omegad);
% Mp = exp(-zeta*omegan*tp)*(cosh(sqrt(zeta^2-1)*omegan*tp)+sinh(sqrt(zeta^2-1)*omegan*tp)*(zeta*omegan-1)/(sqrt(zeta^2-1)*omegan));

y = abs(tr-2)