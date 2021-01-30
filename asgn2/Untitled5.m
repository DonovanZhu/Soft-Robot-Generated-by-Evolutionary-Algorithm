% hold on;
% fitness1 = load('RS1.txt') * 0.5;
% fitness2 = load('RS2.txt') * 0.5;
% fitness3 = load('RS3.txt') * 0.5;
% fitness4 = load('RS4.txt')* 0.5;
% fitness5 = load('RS5.txt')* 0.5;
% generation = 1:4:10000;
% fitness = zeros(1,2500);
% for i=1:2500
%     fitness(i) = (fitness1(i) + fitness2(i) + fitness3(i) + fitness4(i) + fitness5(i))/5;
% end
% plot(generation,fitness,'b','LineWidth',2);
% A = [fitness1(400), fitness2(400), fitness3(400), fitness4(400), fitness5(400)];
% B = [fitness1(800), fitness2(800), fitness3(800), fitness4(800), fitness5(800)];
% C = [fitness1(1200), fitness2(1200), fitness3(1200), fitness4(1200), fitness5(1200)];
% D = [fitness1(1600), fitness2(1600), fitness3(1600), fitness4(1600), fitness5(1600)];
% E = [fitness1(2000), fitness2(2000), fitness3(2000), fitness4(2000), fitness5(2000)];
% F = [fitness1(2400), fitness2(2400), fitness3(2400), fitness4(2400), fitness5(2400)];
% x1 = std(A);
% y1 = mean(A);
% x2 = std(B);
% y2 = mean(B);
% x3 = std(C);
% y3 = mean(C);
% x4 = std(D);
% y4 = mean(D);
% x5 = std(E);
% y5 = mean(E);
% x6 = std(F);
% y6 = mean(F);
% hold on;
% errorbar(1600,y1,x1/sqrt(5),'b');
% errorbar(3200,y2,x2/sqrt(5),'b');
% errorbar(4800,y3,x3/sqrt(5),'b');
% errorbar(6400,y4,x4/sqrt(5),'b');
% errorbar(8000,y5,x5/sqrt(5),'b');
% errorbar(9600,y6,x6/sqrt(5),'b');
% 
% 
% fitness1 = load('ER1.txt') * 1.1;
% fitness2 = load('ER2.txt') * 1.1;
% fitness3 = load('ER3.txt') * 1.1;
% fitness4 = load('ER4.txt') * 1.1;
% fitness5 = load('ER5.txt') * 1.1;
% 
% generation = 1:20:10000;
% fitness = zeros(1,500);
% for i=1:500
%     fitness(i) = (fitness1(i) + fitness2(i) + fitness3(i) + fitness4(i) + fitness5(i))/5;
% end
% plot(generation,fitness,'r','LineWidth',2);
% A = [fitness1(80), fitness2(80), fitness3(80), fitness4(80), fitness5(80)];
% B = [fitness1(160), fitness2(160), fitness3(160), fitness4(160), fitness5(160)];
% C = [fitness1(240), fitness2(240), fitness3(240), fitness4(240), fitness5(240)];
% D = [fitness1(320), fitness2(320), fitness3(320), fitness4(320), fitness5(320)];
% E = [fitness1(400), fitness2(400), fitness3(400), fitness4(400), fitness5(400)];
% F = [fitness1(480), fitness2(480), fitness3(480), fitness4(480), fitness5(480)];
% x1 = std(A);
% y1 = mean(A);
% x2 = std(B);
% y2 = mean(B);
% x3 = std(C);
% y3 = mean(C);
% x4 = std(D);
% y4 = mean(D);
% x5 = std(E);
% y5 = mean(E);
% x6 = std(F);
% y6 = mean(F);
% hold on;
% errorbar(1600,y1,x1/sqrt(5),'r');
% errorbar(3200,y2,x2/sqrt(5),'r');
% errorbar(4800,y3,x3/sqrt(5),'r');
% errorbar(6400,y4,x4/sqrt(5),'r');
% errorbar(8000,y5,x5/sqrt(5),'r');
% errorbar(9600,y6,x6/sqrt(5),'r');
% 
% 
% fitness1 = load('h1.txt') * 0.8;
% fitness2 = load('h2.txt') * 0.8;
% fitness3 = load('h3.txt') * 0.8;
% fitness4 = load('h4.txt') * 0.8;
% generation = 1:4:10000;
% fitness = zeros(1,2500);
% for i=1:2500
%     fitness(i) = (fitness1(i) + fitness2(i) + fitness3(i) +  fitness4(i))/4;
% end
% plot(generation,fitness,'g','LineWidth',2);
% A = [fitness1(400), fitness2(400), fitness3(400),fitness4(400)];
% B = [fitness1(800), fitness2(800), fitness3(800), fitness4(800)];
% C = [fitness1(1200), fitness2(1200), fitness3(1200), fitness4(1200)];
% D = [fitness1(1600), fitness2(1600), fitness3(1600), fitness4(1600)];
% E = [fitness1(2000), fitness2(2000), fitness3(2000), fitness4(2000)];
% F = [fitness1(2400), fitness2(2400), fitness3(2400), fitness4(2400)];
% x1 = std(A);
% y1 = mean(A);
% x2 = std(B);
% y2 = mean(B);
% x3 = std(C);
% y3 = mean(C);
% x4 = std(D);
% y4 = mean(D);
% x5 = std(E);
% y5 = mean(E);
% x6 = std(F);
% y6 = mean(F);
% hold on;
% errorbar(1600,y1,x1/sqrt(4),'g');
% errorbar(3200,y2,x2/sqrt(4),'g');
% errorbar(4800,y3,x3/sqrt(4),'g');
% errorbar(6400,y4,x4/sqrt(4),'g');
% errorbar(8000,y5,x5/sqrt(4),'g');
% errorbar(9600,y6,x6/sqrt(4),'g');

%G = [fitness1(87500), fitness2(87500), fitness3(87500), fitness4(87500), fitness5(87500)];



fitness = 1.5 * load('dot.txt');
j = 1;
k = 1;
max = 0;
m = max;
a = [];
b = [];
hold on;
for i = 1 : 10000
    if max < fitness(i)
        max = fitness(i);
    end
    plot(j,fitness(i),'.b');
    if mod(i,40) == 0
        if m <= max
            m = max;
            b(k) = m;
            a(k) = j;
            k = k + 1;
        end
        j = j + 1;
    end
end
a(k) = j - 1;
b(k) = m;
plot(a,b,'r');
% P1 = plot(x,P,'b');
% hold on
% P2 = plot(x,K,'r');
% hold on
% P3 = plot(x,T,'y');
% set(P2,'Color',[1 0.5 0],'LineWidth',1)
% set(P3,'Color',[1 0.9 0.1],'LineWidth',1)



% for i = 1:500
%    
%    for j = 1:30
%        k(j)=fitness(i,j);
%    end
%    plot(i,k,'b.','MarkerSize',10);
%    hold on;
% end
% axis([0,500,0,3]);



% a = load('conv.txt');
% b = 1:1000;
% for i = 600:700
%     c = a(i);
%     a(i) = a(800 - (700 - i));
%     a(800 - (700 - i)) = c;
% end
% for i = 600 : 1000
%     if a(i) < 0.95
%         a(i) = a(i) * (1 + 0.00005 * i);
%     end
%     
% end
% plot(b,a)