function plotMechanism(q,t)
lDE = 0.15;
lGF = 0.1;
h = 0.5;
Sx1    = q(1);
Sy1    = q(2);
theta1  = q(3);
Sx2    = q(4);
Sy2    = q(5);
theta2 = q(6);
r1 = exp(1i*theta1);
r2 = exp(1i*theta2);

S1 = Sx1+1i*Sy1;
D = S1-lDE/2*r1;
E = S1+lDE/2*r1;

S2 = Sx2+1i*Sy2;
G = S2-lGF/2*r2;
F = S2+lGF/2*r2;
plotComplexPoints(D,E);
hold on
plotComplexPoints(G,F);

C = 1i*h;
plotComplexPoint(C,'k');
plotComplexPoint(D,'k');

hold off
axis equal
axis([-0.25,0.25,0.25,0.75])
grid on
title(['2-bar mechanism simulation, t = ' num2str(t)])
end