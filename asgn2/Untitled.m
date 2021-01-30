fitness1 = load('Random1.txt')/1000;
fitness2 = load('Random2.txt')/1000;
fitness3 = load('Random3.txt')/1000;
fitness4 = load('Random4.txt')/1000;
fitness5 = load('Random5.txt')/1000;
generation = 1:100000;
fitness = zeros(1,100000);
for i=1:100000
    fitness(i) = (fitness1(i) + fitness2(i) + fitness3(i) + fitness4(i) + fitness5(i))/5;
end
A = [fitness1(12500), fitness2(12500), fitness3(12500), fitness4(12500), fitness5(12500)];
B = [fitness1(25000), fitness2(25000), fitness3(25000), fitness4(25000), fitness5(25000)];
C = [fitness1(37500), fitness2(37500), fitness3(37500), fitness4(37500), fitness5(37500)];
D = [fitness1(50000), fitness2(50000), fitness3(50000), fitness4(50000), fitness5(50000)];
E = [fitness1(62500), fitness2(62500), fitness3(62500), fitness4(62500), fitness5(62500)];
F = [fitness1(75000), fitness2(75000), fitness3(75000), fitness4(75000), fitness5(75000)];
G = [fitness1(87500), fitness2(87500), fitness3(87500), fitness4(87500), fitness5(87500)];
x1 = std(A);
y1 = mean(A);
x2 = std(B);
y2 = mean(B);
x3 = std(C);
y3 = mean(C);
x4 = std(D);
y4 = mean(D);
x5 = std(E);
y5 = mean(E);
x6 = std(F);
y6 = mean(F);
x7 = std(G);
y7 = mean(G)
hold on;
errorbar(12500,y1,x1/sqrt(5),'g');
% errorbar(25000,y2,x2/sqrt(5),'g');
% errorbar(37500,y3,x3/sqrt(5),'g');
% errorbar(50000,y4,x4/sqrt(5),'g');
% errorbar(62500,y5,x5/sqrt(5),'g');
% errorbar(75000,y6,x6/sqrt(5),'g');
% errorbar(87500,y7,x7/sqrt(5),'g');
% hold on
% plot(generation,fitness,'g','LineWidth',2);
% axis([0,100000,0,1.5]);
% xlabel('Evaluations');
% ylabel('Fitness');


%  x1 = xlsread('Point.xlsx','A1:A1000');
%  y1 = xlsread('Point.xlsx','B1:B1000');
%  y2 = zeros(1,1000);
%  x2 = zeros(1,1000);
%  a = 10 / 1000;
%  for i = 1 : 1000
%      x2(i) = a * i;
%  end
%  x = 0;
%  for i = 1 : 1000
%     x = x1(i);
%     y2(i) = (((7.3369701223-sin((sin(x)+cos(x))))-sin((sin(x)+cos(x))))*0.000000000);
% %    y2(i) = 2.5*(0.7+sin(x-2.711767314084740))*sin(x);
% %    y2(i) = (2.4983520005*((sin(sin(cos((-1.6981261635-4.1125370036))))+sin((x+(sin(2.9689779351)-2.8835261086))))*sin(x)));
% %    y2(i) = (2.5040284432*((sin(sin(cos((5.8817560351-6.3797875912))))+sin((x+(sin(6.0779442733)-2.5040284432))))*sin(x)));
%     %   y2(i) = ((((-1.7083346050+sin(sin(x)))-sin(((x+x)+sin(-1.7083346050))))-sin(sin(((-0.1948988311+sin(x))+sin(-1.7083346050)))))+sin(x));
%  end
% %  hold on;
% %  plot(x1,y1);
% %hold on;
% %plot(x2,y2,'r','LineWidth',3);
% hold on;
% plot(x1,y1,'r.','MarkerSize',13);
% hold on;
% plot(x1,y2,'b.','MarkerSize',13);
% 
% 
% error = 0;
% for i = 1 : 1000
%    error = error + abs(y1(i) - y2(i)); 
% end
% error