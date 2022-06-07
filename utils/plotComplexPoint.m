function plotComplexPoint(p,color,marker_size)
if nargin < 3
    marker_size = 5; 
    if nargin < 2
        color = 'r';
    end
end
plot(real(p),imag(p),'or',...
    'MarkerFaceColor',color,...
    'MarkerSize',marker_size)
end