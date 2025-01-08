function [q1,q2] = qfunction(wn)
%QFUNCTION Summary of this function goes here
%   Detailed explanation goes here
% Tasarımda ts 2sn olarak seçildiği için ts=4/(zeta*wn) formülünden zeta
% ve wn aralıkları aşağıdaki gibi belirlenir:
% zeta: [0.1 0.99] ve wn: [2.02 20]
% zeta*wn=2 dir.

p=1;
%zeta*wn=2
%q1=p*((((wn^2+2*zeta*wn-1)^2)/4)+wn^2+2*zeta*wn-1);
q1=p*((((wn^2+3)^2)/4)+wn^2+3);

%q2=p*((((-wn^2+2*zeta*wn-1-2*p)^2)/4)-(wn^2+2*zeta*wn-1));
q2=p*((((-wn^2+1)^2)/4)-(wn^2+3));

end