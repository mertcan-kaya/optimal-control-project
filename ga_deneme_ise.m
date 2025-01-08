clc, clear all, close all

A = [];
b = [];
Aeq = [];
beq = [];
% lb = [];
ub = [];

lb = [0.00001 0.00001 0.00001];

% out = ga(@ISEfunc2,3,A,b,Aeq,beq,lb,ub)
% out = ga(@xTQxuTPu,3,A,b,Aeq,beq,lb,ub)
% out = ga(@ITSEfunc,3,A,b,Aeq,beq,lb,ub)
% out = ga(@costTrTpfunc,3,A,b,Aeq,beq,lb,ub)
% out = ga(@costMpfunc,3,A,b,Aeq,beq,lb,ub)
% out = ga(@costTrMpfunc,3,A,b,Aeq,beq,lb,ub)
% out = ga(@costTrMpfunc2,3,A,b,Aeq,beq,lb,ub)
% out = ga(@costTrTpMpfunc,3,A,b,Aeq,beq,lb,ub)
% out = ga(@ISEUfunc,3,A,b,Aeq,beq,lb,ub)
out = ga(@isq1xpufunc,3,A,b,Aeq,beq,lb,ub)

% ISEfunc2(out)

A = [0,1;-1,0];
B = [0;1];

Q = [out(1),0;0,out(2)]
P = out(3)

R = icare(A,B,Q,P,[],[],[]);

K = R*B/P

k1 = K(1);
k2 = K(2);

syms t

x(t) = (exp(-t*(k2/2 + (k2^2 - 4*k1 - 4)^(1/2)/2))*((k2^2 - 4*k1 - 4)^(1/2) - k2 + 2))/(2*(k2^2 - 4*k1 - 4)^(1/2)) + (exp(-t*(k2/2 - (k2^2 - 4*k1 - 4)^(1/2)/2))*(k2 + (k2^2 - 4*k1 - 4)^(1/2) - 2))/(2*(k2^2 - 4*k1 - 4)^(1/2));

Dx(t) = (exp(-t*(k2/2 + (k2^2 - 4*k1 - 4)^(1/2)/2))*(2*k1 - k2 - (k2^2 - 4*k1 - 4)^(1/2) + 2))/(2*(k2^2 - 4*k1 - 4)^(1/2)) - (exp(-(t*(k2 - (k2^2 - 4*k1 - 4)^(1/2)))/2)*(2*k1 - k2 + (k2^2 - 4*k1 - 4)^(1/2) + 2))/(2*(k2^2 - 4*k1 - 4)^(1/2));

u(t) = -K'*[x(t);Dx(t)];

figure
hold on
fplot(u,[0 10])
fplot(x,[0 10])
fplot(Dx,[0 10])
xlabel('Time (s)')
legend('u','x_1','x_2')

figure
fplot(-x,[0 10])
xlabel('Time (s)')
ylabel('Error')