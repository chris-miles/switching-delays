clear all
close all


figure('Position',[0,0,500, 200]);


rng(1);


I = 10;
K = 9.5;
d = 1;
nn = 4;
M = 2;

epsilon = 10;


Tmax = 50;

nTau = 2;
%colormap turbo
addpath('cbrewer/')
colmap = cbrewer('seq','GnBu',3);
cc = colmap(end-nTau:end,:);

tauvals = fliplr(linspace(0.6, 1.2, nTau));
for n = 1:nTau
    tau = tauvals(n);
    ccc = cc(n,:);

p0 = 0; % initial value
history0 = @(t) p0*ones(1,1); % pick (constant) history
effective_DDE = @(t,y,z)  (I-d*y(1) - K*(z(1,1)^nn)/(M^nn+z(1,1)^nn));
sol_1delay = dde23(effective_DDE, tau, history0, [0, Tmax]);

plot(sol_1delay.x,sol_1delay.y,'LineWidth',2.5,'color',ccc)
hold on;
end 

xlim([-1, Tmax]);
ylim([0 8]);
xlabel('t')
ylabel('y(t)');
set(gca,'LineWidth',1.5)
set(gca,'FontSize',12)
box off;
%pbaspect([4 3 1])

