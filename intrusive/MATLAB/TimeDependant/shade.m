

x = linspace(0,1);
y1 = 1 + exp(-4*x);
y2 = 1 - exp(-4*x);
figure(1)
plot(x, y1)
hold on
plot(x, y2)
% patch([x x], [y1 fliplr(y2)], 'b')
xn = [x fliplr(x)];
yn = [y1, fliplr(y2)];
fill(xn, yn, 'b');
hold off