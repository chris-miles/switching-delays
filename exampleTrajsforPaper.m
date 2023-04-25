clear all
close all


rng(1);

I = 10;
K = 9.5;
d = 1;
nn = 4;
M = 2;

epsilon = 10;


Tmax = 50;
tau1=1;
tau2=3;
alpha=1;
beta=1;

[y,t,jump_times,states] = run_1stochsim(tau1,tau2,alpha,beta,epsilon,I, K, d, nn, M,Tmax);

jump_1to2 = jump_times(1:2:end);
jump_2to1 = jump_times(2:2:end);


figure('Position',[0,0,500, 200]);
plot(t,y,'LineWidth',2);
hold on;
xline(jump_1to2,'red','LineWidth',1.5,'LineStyle','--');
xline(jump_2to1,'blue','LineWidth',1.5,'LineStyle','--');


xlim([-1, Tmax]);
xlabel('t')
ylabel('y(t)');
set(gca,'LineWidth',1.5)
set(gca,'FontSize',15)
box off;



%% second figure showing eps dependence 
nsims = 3;


figure('position',[0,0, 500, 500]);

subplot(4,1,1);
epsilon = 5;
for n=1:nsims
[y,t,jump_times,states] = run_1stochsim(tau1,tau2,alpha,beta,epsilon,I, K, d, nn, M,Tmax);
plot(t,y,'LineWidth',2);
hold on;
end 

xlim([-1, Tmax]);
xlabel('t')
ylabel('y(t)');
set(gca,'LineWidth',1.5)
set(gca,'FontSize',12)
box off;
title('$\varepsilon=5$','interpreter','LaTeX')


subplot(4,1,2);
epsilon = 1;
for n=1:nsims
[y,t,jump_times,states] = run_1stochsim(tau1,tau2,alpha,beta,epsilon,I, K, d, nn, M,Tmax);
plot(t,y,'LineWidth',2);
hold on;
end 

xlim([-1, Tmax]);
xlabel('t')
ylabel('y(t)');
set(gca,'LineWidth',1.5)
set(gca,'FontSize',12)
box off;
title('$\varepsilon=1$','interpreter','LaTeX')



subplot(4,1,3);
epsilon = .1;
for n=1:nsims
[y,t,jump_times,states] = run_1stochsim(tau1,tau2,alpha,beta,epsilon,I, K, d, nn, M,Tmax);
plot(t,y,'LineWidth',2);
hold on;
end 

xlim([-1, Tmax]);
xlabel('t')
ylabel('y(t)');
set(gca,'LineWidth',1.5)
set(gca,'FontSize',12)
box off;
title('$\varepsilon=.1$','interpreter','LaTeX')
ylim([0 10]);



subplot(4,1,4);
epsilon = .05;
for n=1:nsims
[y,t,jump_times,states] = run_1stochsim(tau1,tau2,alpha,beta,epsilon,I, K, d, nn, M,Tmax);
plot(t,y,'LineWidth',2);
hold on;
end 


p1 = alpha/(alpha+beta);
p2 = 1-p1;


p0 = 0; % initial value
history0 = @(t) p0*ones(1,1); % pick (constant) history
effective_DDE = @(t,y,z)  p1*(I-d*y(1) - K*(z(1,1)^nn)/(M^nn+z(1,1)^nn))+...
    p2*(I-d*y(1) - K*(z(1,2)^nn)/(M^nn+z(1,2)^nn));
eff_sol = dde23(effective_DDE, [tau1, tau2], history0, [0, Tmax]);

plot(eff_sol.x,eff_sol.y,'black','LineStyle','--','LineWidth',1.5)


xlim([-1, Tmax]);
ylim([0 10]);
xlabel('t')
ylabel('y(t)');
set(gca,'LineWidth',1.5)
set(gca,'FontSize',12)
box off;
title('$\varepsilon=0.01$','interpreter','LaTeX')

