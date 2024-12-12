function handle = stairs(ax, x, y, varargin)

x_new = zeros(1, length(x)*2-1);
y_new = zeros(1, length(y)*2-1);

for k = 1:length(x)
    if k == length(x)
        x_new(2*k-1) = x(k);
        y_new(2*k-1) = y(k);
        break
    end

    x_new(2*k-1) = x(k);
    x_new(2*k) = x(k+1);
    y_new(2*k-1) = y(k);
    y_new(2*k) = y(k);
end

handle = plot(ax, x_new, y_new, varargin{:});
end