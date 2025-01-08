clc, clear all, close all

syms s w t

Gs = 1/(s^2+1);
K = [0.4142 1.3522];

eqn = s^2+K(2)*s+K(1)+1;

eqt = ilaplace(1/eqn)

dt = 0.1; time = 0:dt:10;

num = zeros(1,length(time));
for i = 1:length(time)
    num(:,i) = vpa(subs(eqt,t,time(i)));
end

figure
plot(time,num)
% Xs = [x1s;x2s];
% 
% Us = -K*Xs;
% Ys = Gs*Us;
% Yt = ilaplace(Ys);
% 
% dt = 0.1; time = 0:dt:120;
% 
% unum = zeros(2,length(time));
% ynum = zeros(2,length(time));
% for i = 2:length(time)
%     unum(:,i) = vpa(subs(ut1,t,time(i)));
%     ynum(:,i) = vpa(subs(yt1,t,time(i)));
% end
% 
% figure
% subplot(2,1,1)
% plot(time,unum1)
% xlabel('time (s)')
% ylabel('u_1(t)')
% subplot(2,1,2)
% plot(time,unum2)
% xlabel('time (s)')
% ylabel('u_2(t)')
