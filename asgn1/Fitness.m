x = xlsread('Point0.xlsx','A1:A1000');
y = xlsread('Point0.xlsx','B1:B1000');
fitness = 0;
for i = 1:999
    xd = x(i+1) - x(i);
    yd = y(i+1) - y(i);
    fitness = fitness + sqrt(xd^2 + yd^2);
end