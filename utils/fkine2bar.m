function S = fkine2bar(theta,l1,l2)
lDE = 0.15;
lGF = 0.1;
h = 0.5;
theta1 = theta(1);
theta2 = theta(2);
S(1,1) = (0.5*lDE-l1)*cos(theta1);
S(1,2) = (0.5*lDE-l1)*sin(theta1)+h;
S(1,3) = (0.5*lGF-l2)*cos(theta2)+S(1,1)-0.5*lDE*cos(theta1);
S(1,4) = (0.5*lGF-l2)*sin(theta2)+S(1,2)-0.5*lDE*sin(theta1);
end