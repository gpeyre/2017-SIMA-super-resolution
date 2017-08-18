function plot_measure_2d(x,a,col,scaling)

if nargin<4
    scaling = 50;
end
if nargin<3
    col = 'b';
end

x = x(:);
for i=1:length(a)
    plot(x(i), x(end/2+i), '.', 'color', col, 'MarkerSize', scaling*max(1e-2,a(i)));
end

end