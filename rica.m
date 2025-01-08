clc, clear all

A = [0,1;-1,0];
B = [0;1];

Q = [1,0;0,1];
P = 1;

R = icare(A,B,Q,P,[],[],[]);

K = R*B;

k1 = K(1);
k2 = K(2);

model = 'opti_ext';

open_system(model)
load_system(model)
out = sim(model,'StartTime','0','StopTime','10');

figure
plot(out.tout,out.u.Data,out.tout,out.x1.Data,out.tout,out.x2.Data)
xlabel('Time (s)')
legend('u','x1','x2')