clc, clear all

%% Genetic Algoritm

A = [];
b = [];
Aeq = [];
beq = [];
% lb = [];
ub = [];

lb = [0.00001 0.00001 0.00001];

% Select a cost function
% out = ga(@ISEfunc,3,A,b,Aeq,beq,lb,ub) % ISE: integral e^2 dt t:0,10 or inf
% out = ga(@xTQxuTPu,3,A,b,Aeq,beq,lb,ub) % LQ: integral x'Qx+u'P2 dt t:0,10 or inf
% out = ga(@ITSEfunc,3,A,b,Aeq,beq,lb,ub) % ITSE: integral t*e^2 dt t:0,10 or inf
% out = ga(@ISEUfunc,3,A,b,Aeq,beq,lb,ub) % ISEU: integral e^2+u^2 dt t:0,10 or inf
% out = ga(@isq1xpufunc,3,A,b,Aeq,beq,lb,ub) % Cost: integral A*e^2+B*u^2 dt t:0,10
% out = ga(@costTr2func,3,A,b,Aeq,beq,lb,ub) % Cost: tr = 2
% out = ga(@costTr2Mp0func,3,A,b,Aeq,beq,lb,ub) % Cost: tr = 2, Mp = 0
% out = ga(@costTr2_100Mp0func,3,A,b,Aeq,beq,lb,ub) % Cost: tr = 2, Mp = 0 iwth 100x importance

% out = ga(@costTr1func,3,A,b,Aeq,beq,lb,ub) % Cost: tr = 2
% out = ga(@costTr2Mp0func,3,A,b,Aeq,beq,lb,ub) % Cost: tr = 2, Mp = 0
out = ga(@costTr1_100Mp0func,3,A,b,Aeq,beq,lb,ub) % Cost: tr = 2, Mp = 0 iwth 100x importance


%% State Space Representation of the System

A = [0,1;-1,0];
B = [0;1];

%% The Matrix Ricatti Equation

Q = [out(1),0;0,out(2)]
P = out(3)

R = icare(A,B,Q,P,[],[],[]);

K = R*B/P

k1 = K(1);
k2 = K(2);

%% Control Signal and System States

syms t

x(t) = (exp(-t*(k2/2 + (k2^2 - 4*k1 - 4)^(1/2)/2))*((k2^2 - 4*k1 - 4)^(1/2) - k2 + 2))/(2*(k2^2 - 4*k1 - 4)^(1/2)) + (exp(-t*(k2/2 - (k2^2 - 4*k1 - 4)^(1/2)/2))*(k2 + (k2^2 - 4*k1 - 4)^(1/2) - 2))/(2*(k2^2 - 4*k1 - 4)^(1/2));

Dx(t) = (exp(-t*(k2/2 + (k2^2 - 4*k1 - 4)^(1/2)/2))*(2*k1 - k2 - (k2^2 - 4*k1 - 4)^(1/2) + 2))/(2*(k2^2 - 4*k1 - 4)^(1/2)) - (exp(-(t*(k2 - (k2^2 - 4*k1 - 4)^(1/2)))/2)*(2*k1 - k2 + (k2^2 - 4*k1 - 4)^(1/2) + 2))/(2*(k2^2 - 4*k1 - 4)^(1/2));

u(t) = -K'*[x(t);Dx(t)];

% Natural frequency
omegan = sqrt((k1 + 1))
% Damping ratio
zeta = k2/(2*omegan)

% Damped natural frequency
omegad = omegan*sqrt(1-zeta^2);
% Phase shift
phi = acos(zeta);

% Time shift caused by initial conditions
tshft = -acosh((omegan-zeta)/sqrt(1-2*zeta*omegan+omegan^2))/(sqrt(zeta^2-1)*omegan);

% Rise time
tr = tshft+real((pi-phi)/omegad)
% Rise point
Mr = -exp(-zeta*omegan*tr)*(cosh(sqrt(zeta^2-1)*omegan*tr)+sinh(sqrt(zeta^2-1)*omegan*tr)*(zeta*omegan-1)/(sqrt(zeta^2-1)*omegan));
% Peak time
tp = tshft+real(pi/omegad);
% Peak overshoot
Mp = -exp(-zeta*omegan*tp)*(cosh(sqrt(zeta^2-1)*omegan*tp)+sinh(sqrt(zeta^2-1)*omegan*tp)*(zeta*omegan-1)/(sqrt(zeta^2-1)*omegan))

q1 = Q(1,1);
q2 = Q(2,2);
p = P;

cost10 = double(int(q1*x*x+q2*Dx*Dx+p*u*u,t,0,10))
costInf = double(int(q1*x*x+q2*Dx*Dx+p*u*u,t,0,inf))

figure
hold on
fplot(u,[0 10])
fplot(x,[0 10])
fplot(Dx,[0 10])
xlabel('Time [s]')
legend('u','x_1','x_2')

figure
hold on
fplot(-x,[0 10])
plot(tr,Mr,'o',tp,Mp,'o')
yline(0,':')
xline(tr,':')
yline(Mp,':')
xline(tp,':')
xlabel('Time [s]')
ylabel('e(t)')
legend('Error','Rise','Peak','Location','southeast')