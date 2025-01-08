clc, clear all, close all

syms x(t) k1 k2

eqn = diff(x,t,2) + k2*diff(x,t) + (k1 + 1)*x == 0;

Dx = diff(x,t);
cond = [x(0)==1, Dx(0)==-1];

xSol(t) = dsolve(eqn,cond);
DxSol(t) = diff(xSol,t);

u(t) = -k1*xSol(t) - k2*DxSol(t);
% e = -xSol(t);
% 
% ISE = double(int(e^2,t,[0 10]))
% ITSE = double(int(t*e^2,t,[0 10]))
% 
% % solve_e = double(solve(e,t));
% % tr = real(solve_e(solve_e>0))
% tr = double(vpasolve(e))
% Ptr = double(subs(e,t,tr));
% 
% % solve_dx = double(solve(DxSol,t));
% % tMp = real(solve_dx(solve_dx>0))
% tMp = double(vpasolve(DxSol,[0 10]))
% Mp = double(subs(e,t,tMp))
% 
% pct5 = 0.005;
% ts5 = double(vpasolve(abs(e)==pct5,[tMp 10]))
% Pts5 = double(subs(e,t,ts5));
% pct2 = 0.002;
% ts2 = double(vpasolve(abs(e)==pct2,[tMp 10]))
% Pts2 = double(subs(e,t,ts5));
% 
% figure
% fplot(u,[0 10]);
% hold on
% fplot(xSol,[0 10]);
% fplot(DxSol,[0 10]);
% axis([0 10 -1.5 1.5])
% xlabel('Time (s)')
% legend('u','x_1','x_2')
% 
% figure
% fplot(e,[0 10]);
% hold on
% plot(tr,Ptr,'o',tMp,Mp,'o',ts5,Pts5,'o')
% yline(0,':')
% xline(tr,':')
% yline(Mp,':')
% xline(tMp,':')
% yline(Pts5,':')
% xline(ts5,':')
% axis([0 10 -1 0.25])
% xlabel('Time (s)')
% ylabel('e(t)')
% legend('Error','Rise','Peak','Settling','Location','southeast')