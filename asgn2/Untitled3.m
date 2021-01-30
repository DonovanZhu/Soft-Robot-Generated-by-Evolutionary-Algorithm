fitness1 = load('fitness1_GP.txt');
fitness2 = load('fitness2_GP.txt');
fitness3 = load('fitness3_GP.txt');
fitness4 = load('fitness4_GP.txt');
fitness5 = load('fitness5_GP.txt');
generation = 1:100:100000;
fitness = zeros(1,1000);
for i=1:1000
    fitness(i) = (fitness1(i) + fitness2(i) + fitness3(i) + fitness4(i) + fitness5(i))/5;
end
A = [fitness1(125), fitness2(125), fitness3(125), fitness4(125), fitness5(125)];
B = [fitness1(250), fitness2(250), fitness3(250), fitness4(250), fitness5(250)];
C = [fitness1(375), fitness2(375), fitness3(375), fitness4(375), fitness5(375)];
D = [fitness1(500), fitness2(500), fitness3(500), fitness4(500), fitness5(500)];
E = [fitness1(625), fitness2(625), fitness3(625), fitness4(625), fitness5(625)];
F = [fitness1(750), fitness2(750), fitness3(750), fitness4(750), fitness5(750)];
G = [fitness1(875), fitness2(875), fitness3(875), fitness4(875), fitness5(875)];
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
errorbar(12500,y1,x1/sqrt(5),'c');
% errorbar(25000,y2,x2/sqrt(5),'c');
% errorbar(37500,y3,x3/sqrt(5),'c');
% errorbar(50000,y4,x4/sqrt(5),'c');
% errorbar(62500,y5,x5/sqrt(5),'c');
% errorbar(75000,y6,x6/sqrt(5),'c');
% errorbar(87500,y7,x7/sqrt(5),'c');
% hold on
% plot(generation,fitness,'c','LineWidth',2);
% axis([0,100000,0,1.5]);
% xlabel('Evaluations');
% ylabel('Fitness');