function plot2dPoints(p1,p2,dotted,color)
if nargin < 4
    color = 'k';
    if nargin < 3
        dotted = 0;
    end
end
x = [p1(1),p2(1)];
y = [p1(2),p2(2)];
if dotted
    plot(x,y,['--' color]);
else
    plot(x,y,color);
end
axis equal
end

