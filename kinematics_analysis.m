clear all;close all;
%% set geometric parameters
syms lDE lGF l1 l2 h
%% position constraints
% general coordinates
% number of movable link
syms Sx1 Sy1 theta1 Sx2 Sy2 theta2
q(1,1) = Sx1;
q(2,1) = Sy1;
q(3,1) = theta1;
q(4,1) = Sx2;
q(5,1) = Sy2;
q(6,1) = theta2;
n = size(q,1);

% contraint loops
% twice number of constraints
Phi(1,1) = Sx1-(0.5*lDE-l1)*cos(theta1);
Phi(2,1) = Sy1-(0.5*lDE-l1)*sin(theta1)-h;
Phi(3,1) = Sx2-(0.5*lGF-l2)*cos(theta2)-Sx1+0.5*lDE*cos(theta1);
Phi(4,1) = Sy2-(0.5*lGF-l2)*sin(theta2)-Sy1+0.5*lDE*sin(theta1);
s = size(Phi,1);

%% velocity constraints
for i = 1:n
    for j = 1:s
        Phi_q(j,i) = diff(Phi(j),q(i),1);
    end
end

%% acceleration constraints
syms dSx1 dSy1 dtheta1 dSx2 dSy2 dtheta2
dq(1,1) = dSx1;
dq(2,1) = dSy1;
dq(3,1) = dtheta1;
dq(4,1) = dSx2;
dq(5,1) = dSy2;
dq(6,1) = dtheta2;
tmpv = Phi_q*dq;
for i = 1:n
    for j = 1:s
        tmp(j,i) = diff(tmpv(j),q(i),1);
    end
end
gamma = -tmp*dq;