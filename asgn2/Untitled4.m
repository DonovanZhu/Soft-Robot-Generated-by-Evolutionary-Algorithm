 x1 = xlsread('Point.xlsx','A1:A1000');
 y1 = xlsread('Point.xlsx','B1:B1000');
 y2 = zeros(1,1000);
 x2 = zeros(1,1000);
 a = 10 / 1000;
 for i = 1 : 1000
     x2(i) = a * i;
 end
 x = 0;
 for i = 1 : 1000
     x = x2(i);
     y2(i) = 2.5*(0.7+sin(x-2.711767314084740))*sin(x);
     y2(i) = 2.0187994018*sin(x) + 1.4767801706;
%    y2(i) = (2.5040284432*((sin(sin(cos((5.8817560351-6.3797875912))))+sin((x+(sin(6.0779442733)-2.5040284432))))*sin(x)));
    %   y2(i) = ((((-1.7083346050+sin(sin(x)))-sin(((x+x)+sin(-1.7083346050))))-sin(sin(((-0.1948988311+sin(x))+sin(-1.7083346050)))))+sin(x));
 end
%  hold on;
%  plot(x1,y1);

hold on;
plot(x1,y1,'b.','MarkerSize',10);
hold on;
plot(x2,y2,'r','LineWidth',1);
% hold on;
% plot(x1,y2,'b.','MarkerSize',13);


error = 0;
for i = 1 : 1000
   error = error + abs(y1(i) - y2(i)); 
end
error;