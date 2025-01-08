clc, clear all, close all

syms x(t)

A = [0,1;-1,0];
B = [0;1];

Q = [1,0;0,1];
P = 1;

% R = icare(A,B,Q,P,[],[],[]);
R = care(A,B,Q,P,[],[]);

K = R*B;

eqn = diff(x,t,2) + K(2)*diff(x,t) + (K(1) + 1)*x == 0;

omegan = sqrt((K(1) + 1))
zeta = K(2)/(2*omegan)

omegad = omegan*sqrt(1-zeta^2)
theta = acos(zeta)

trise = (pi-theta)/omegad

tpeak = pi/omegad

tsettling5 = 3/(zeta*omegan)
tsettling2 = 4/(zeta*omegan)

% KT = [-1+omegan^2,2*zeta*omegan]

Dx = diff(x,t);
cond = [x(0)==1, Dx(0)==-1];

xSol(t) = dsolve(eqn,cond);
DxSol(t) = diff(xSol,t);

u(t) = -K(1)*xSol(t) - K(2)*DxSol(t);
e = -xSol(t);

% isefunc = int(e^2,t,[0 10])
% ISE = double(int(e^2,t,[0 10]))
% ITSE = double(int(t*e^2,t,[0 10]))

what = double(solve(xSol==-DxSol,t))

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

trise2 = tshift+(pi-theta)/omegad
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

expoP = +exp(eigreal.*time)/sqrt(1-zeta^2);
expoN = -exp(eigreal.*time)/sqrt(1-zeta^2);

figure
fplot(u,[-2 10]);
hold on
fplot(xSol,[-2 10]);
fplot(DxSol,[-2 10]);
plot(time,expoP,time,expoN)
axis([-2 10 -2 2])
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