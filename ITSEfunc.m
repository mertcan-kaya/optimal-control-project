function y = ITSEfunc(G)

A = [0,1;-1,0];
B = [0;1];

Q = [G(1),0;0,G(2)];
P = G(3);

G
R = icare(A,B,Q,P,[],[],[])

K = R*B/P;

k1 = K(1);
k2 = K(2);

% % a) integrate from 0 to T

% T = 10;

% y = (2*k1^2 - 2*k1*k2 + 6*k1 + k2^4 - 2*k2^3 + k2^2 - 2*k2 + 4)/(4*k2^2*(k1 + 1)^2) - ((4*exp(-T*k2)*(k1 - k2 + 2))/k2^2 + (exp(-T*(k2 - (k2^2 - 4*k1 - 4)^(1/2)))*(2*k1 + 2*k2 - (k2 - 2)*(k2^2 - 4*k1 - 4)^(1/2) - k2^2))/(k2 - (k2^2 - 4*k1 - 4)^(1/2))^2 + T*exp(-T*(k2 - (k2^2 - 4*k1 - 4)^(1/2)))*(k2 + (- 2*k2^2 + 4*k2 + 2*k1)/(k2 - (k2^2 - 4*k1 - 4)^(1/2)) - 2) + (exp(-T*(k2 + (k2^2 - 4*k1 - 4)^(1/2)))*(2*k1 + 2*k2 + (k2 - 2)*(k2^2 - 4*k1 - 4)^(1/2) - k2^2))/(k2 + (k2^2 - 4*k1 - 4)^(1/2))^2 + T*exp(-T*(k2 + (k2^2 - 4*k1 - 4)^(1/2)))*(k2 + (- 2*k2^2 + 4*k2 + 2*k1)/(k2 + (k2^2 - 4*k1 - 4)^(1/2)) - 2) + (4*T*exp(-T*k2)*(k1 - k2 + 2))/k2)/(2*(- k2^2 + 4*k1 + 4))

% b) integrate from 0 to infinite
syms t

x(t) = (exp(-t*(k2/2 + (k2^2 - 4*k1 - 4)^(1/2)/2))*((k2^2 - 4*k1 - 4)^(1/2) - k2 + 2))/(2*(k2^2 - 4*k1 - 4)^(1/2)) + (exp(-t*(k2/2 - (k2^2 - 4*k1 - 4)^(1/2)/2))*(k2 + (k2^2 - 4*k1 - 4)^(1/2) - 2))/(2*(k2^2 - 4*k1 - 4)^(1/2));

y = double(int(t*x^2,t,0,inf))
