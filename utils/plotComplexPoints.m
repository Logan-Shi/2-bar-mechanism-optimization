function plotComplexPoints(p1,p2,dotted,color)
if nargin < 4
    color = 'k';
    if nargin < 3
        dotted = 0;
    end
end
x = [real(p1),real(p2)];
y = [imag(p1),imag(p2)];
if dotted
    plot(x,y,['--' color]);
else
    plot(x,y,color);
end
axis equal
end

