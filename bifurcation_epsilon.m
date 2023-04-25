clear all
close all


I = 10;
K = 9.5;
d = 1;
nn = 4;
M = 2;


% probability of being in each state

alpha=1;
beta=1;

p1 = beta/(alpha+beta);
p2 = 1-p1;


p0 = 0; % initial value
history0 = @(t) p0*ones(1,1); % pick (constant) history
dydt = @(t,y,z)  I-d*y(1) - K*(z(1,1)^nn)/(M^nn+z(1,1)^nn);

effective_DDE = @(t,y,z)  p1*(I-d*y(1) - K*(z(1,1)^nn)/(M^nn+z(1,1)^nn))+...
    p2*(I-d*y(1) - K*(z(1,2)^nn)/(M^nn+z(1,2)^nn));


tau1=1;%1
tau2=3;%3;

nEps = 30;
epsvals = logspace(-1.5,1.5,nEps);
%epsvals = logspace(-1,2,nEps);
nSims = 150; % was 150 in good version

finalVals = zeros(nEps,nSims);
maxVals = zeros(nEps,nSims);
minVals = zeros(nEps,nSims);
fftVals = zeros(nEps,nSims);



%Tmax=500;

%eff_sol = dde23(effective_DDE, [tau1, tau2], history0, [0, Tmax]);
tic;

fracend = .01; %what % to take end 
%tvals_end =  linspace(Tmax*(1-fracend), Tmax,1000);

Tmaxvals =10*epsvals+200;
% or maybe 100*epsvals+100;

for n = 1:nEps
    epsilon = epsvals(n);
    Tmaxn = Tmaxvals(n);
    tvals_end =  linspace(Tmaxn*(1-fracend), Tmaxn,1000);
    
    parfor m = 1:nSims
        [y,t,~,~]= run_1stochsim(tau1,tau2,alpha,beta,epsilon,I, K, d, nn, M,Tmaxn);
        finalVals(n,m) = y(end);
        yy = interp1(t,y,tvals_end);
        
        ttt = t(round(end/2):end);
        yyy = y(round(end/2):end);

        fftval = max(2*abs(nufft(yyy-mean(yyy),ttt))/length(yyy));

        fftVals(n,m) = fftval;
        maxVals(n,m) = max(yy);
        minVals(n,m) = min(yy);

        disp([n,m]);
    end
end


tt = toc;
disp(tt);

%save('big_eps_sweep5')

figure('position',[0,0, 400, 300]);
%scatter(log10(epsvals),finalVals,25,'filled', ...
%    'MarkerFacealpha',17.5/nSims,'MarkerEdgeColor','none','MarkerFaceColor','black');


scatter(log10(epsvals),maxVals,125,'filled', ...
'MarkerFaceAlpha',min([0.5,5/nSims]),'MarkerEdgeColor','none','MarkerFaceColor','black');

hold on;
scatter(log10(epsvals),minVals,125,'filled', ...
'MarkerFaceAlpha',min([0.5,5/nSims]),'MarkerEdgeColor','none','MarkerFaceColor','black');

tau1_sol = dde23(dydt, tau1, history0, [0, Tmaxvals(end)]);
tau2_sol = dde23(dydt, tau2, history0, [0, Tmaxvals(end)]);

max2 = max(tau2_sol.y(round(end/10):end));
min2 = min(tau2_sol.y(round(end/10):end));

max1 = max(tau1_sol.y(round(end/10):end));
min1 = min(tau1_sol.y(round(end/10):end));

yline(min1,'red')
yline(max1,'red')

yline(min2,'blue')
yline(max2,'blue')

pbaspect([ 4 3 1 ]);
ylim([0 10]);

xlabel('$\log_{10}(\varepsilon)$','Interpreter','LaTeX')
ylabel('max/min','Interpreter','LaTeX');
set(gca,'LineWidth',1.5)
set(gca,'FontSize',15)
set(gcf,'Renderer','Painters')

        fftval1 = max(2*abs(nufft(tau1_sol.y-mean(tau1_sol.y),tau1_sol.x))/length(tau1_sol.y));
        fftval2 = max(2*abs(nufft(tau2_sol.y-mean(tau2_sol.y),tau2_sol.x))/length(tau2_sol.y));



figure('position',[0,0, 400, 300]);
scatter(log10(epsvals),fftVals,175,'filled', ...
'MarkerFaceAlpha',min([0.5,5/nSims]),'MarkerEdgeColor','none','MarkerFaceColor','black');
hold on;
%yline(fftval1)
%yline(fftval2)
pbaspect([ 4 3 1 ]);

xlabel('$\log_{10}(\varepsilon)$','Interpreter','LaTeX')
ylabel('\maxS(f)','Interpreter','LaTeX');
set(gca,'LineWidth',1.5)
set(gca,'FontSize',15)
set(gcf,'Renderer','Painters')
