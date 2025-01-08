clc, clear all, close all

syms x(t) T

A = [0,1;-1,0];
B = [0;1];

Q = [1,0;0,1];
P = 1;

R = icare(A,B,Q,P,[],[],[]);

K = R*B/P;

eqn = diff(x,t,2) + K(2)*diff(x,t) + (K(1) + 1)*x == 0;

omegan = sqrt((K(1) + 1));
zeta = K(2)/(2*omegan);

omegad = omegan*sqrt(1-zeta^2);
phi = acos(zeta);

% c(t) = -(exp(-zeta*omegan*t)/sqrt(1-zeta^2))*sin(omegad*t+phi);
x1(t) = exp(-zeta*omegan*t)*(cosh(sqrt(zeta^2-1)*omegan*t)+sinh(sqrt(zeta^2-1)*omegan*t)*(zeta*omegan-1)/(sqrt(zeta^2-1)*omegan));
x2(t) = exp(-zeta*omegan*t)*(-cosh(sqrt(zeta^2-1)*omegan*t)+sinh(sqrt(zeta^2-1)*omegan*t)*(zeta-omegan)/(sqrt(zeta^2-1)));

XX(t) = [x1(t);x2(t)];
u(t) = -K'*XX(t);

figure
fplot(u,[0 10])
hold on
fplot(x1,[0 10])
fplot(x2,[0 10])
axis([0 10 -1.25 1.25])

tshft = -acosh((omegan-zeta)/sqrt(1-2*zeta*omegan+omegan^2))/(sqrt(zeta^2-1)*omegan)

% trise = (pi-theta)/omegad
% 
% tpeak = pi/omegad
% 
% tsettling5 = 3/(zeta*omegan)
% tsettling2 = 4/(zeta*omegan)

Dx = diff(x,t);
cond = [x(0)==1, Dx(0)==-1];

xSol(t) = dsolve(eqn,cond);
DxSol(t) = diff(xSol,t);
xVec(t) = [xSol(t);DxSol(t)];

u(t) = -K'*xVec(t);
e = -xSol(t);

% tf = 10;
% isefunc = int(e^2,t,[0 tf])
% ISE = double(int(e^2,t,[0 tf]))
% ITSE = double(int(t*e^2,t,[0 tf]))

solve_e = double(solve(e,t));
tr = real(solve_e(solve_e>0))
% tr = double(vpasolve(e))
Ptr = double(subs(e,t,tr));

solve_dx = double(solve(DxSol,t))
tshift = real(solve_dx(solve_dx<0))
tMp = real(solve_dx(solve_dx>0))
tMp = double(vpasolve(DxSol,[0 10]))
Mp = double(subs(e,t,tMp))

xshift = double(xSol(tshift))

trise2 = tshift+(pi-phi)/omegad
tpeak2 = tshift+pi/omegad
tsettling52 = tshift-(log(0.05/xshift)/(zeta*omegan))
tsettling52 = tshift+(3/(zeta*omegan))
tsettling22 = tshift-(log(0.02/xshift)/(zeta*omegan))
tsettling22 = tshift+(4/(zeta*omegan))

xsettling5 = double(xSol(tsettling52))
xsettling2 = double(xSol(tsettling22))

pct5 = 0.05;
ts5 = double(vpasolve(abs(e)==pct5,[tMp 10]))
Pts5 = double(subs(e,t,ts5));
pct2 = 0.02;
ts2 = double(vpasolve(abs(e)==pct2,[tMp 10]))
Pts2 = double(subs(e,t,ts5));

time = -2:0.1:10;
eigreal = -zeta*omegan
eigimP = omegan*sqrt(1-zeta^2)

eig1 = complex(eigreal,eigimP)
eig2 = complex(eigreal,-eigimP)

figure
fplot(u,[0 10]);
hold on
fplot(xSol,[0 10]);
fplot(DxSol,[0 10]);
axis([0 10 -1.25 1.25])
xlabel('Time (s)')
legend('u','x_1','x_2')

figure
fplot(e,[0 10]);
hold on
plot(tr,Ptr,'o',tMp,Mp,'o',ts5,Pts5,'o')
yline(0,':')
xline(tr,':')
yline(Mp,':')
xline(tMp,':')
yline(Pts5,':')
xline(ts5,':')
axis([0 10 -1 0.25])
xlabel('Time (s)')
ylabel('e(t)')
legend('Error','Rise','Peak','Settling','Location','southeast')