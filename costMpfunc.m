function y = costMpfunc(G)

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

tshft = -acosh((omegan-zeta)/sqrt(1-2*zeta*omegan+omegan^2))/(sqrt(zeta^2-1)*omegan);

% trise = tshft+real((pi-phi)/omegad)
t = tshft+real(pi/omegad);
Mp = exp(-zeta*omegan*t)*(cosh(sqrt(zeta^2-1)*omegan*t)+sinh(sqrt(zeta^2-1)*omegan*t)*(zeta*omegan-1)/(sqrt(zeta^2-1)*omegan));

y = abs(Mp);