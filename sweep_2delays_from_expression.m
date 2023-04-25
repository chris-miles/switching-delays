clear all
close all


I = 10;
K = 9.5;
d = 1;
nn = 4; 
M = 2;



% probability of being in each state
p1 = 0.5;
p2 = 1-p1;

effect = @(y)  (I-d*y - K*(y^nn)/(M^nn+y^nn)).^2;

[ybar,~,searchflag] = fminsearch(effect,1);


fprime = @(y) ((M^(nn))*nn*y.^(nn-1))./((M^nn+y^nn).^2);

eq1 = @(w,t1,t2) w-K*p1*fprime(ybar).*sin(w*t1)-K*p2*fprime(ybar)*sin(w*t2);
eq2 = @(w,t1,t2) d+K*p1*fprime(ybar).*cos(w*t1)+K*p2*fprime(ybar)*cos(w*t2);


eqtot = @(w,t1,t2) log(eq1(w,t1,t2).^2+eq2(w,t1,t2).^2);


Nsweep = 100;
taumax = 5;
taumin = 0.5;

tauvals = linspace(taumin,taumax,Nsweep);

%t1fix=1;




%[wt2_opt,eqval,searchflag] = fminsearch(eqtot_t1fix,[1,2]);
%wopt = wt2_opt(1)
%t2opt = exp(wt2_opt(2))

%options = optimoptions('simulannealbnd','PlotFcns',...
%         {@saplotbestx,@saplotbestf,@saplotx,@saplotf});

%[wt2_opt,eqval,searchflag] = simulannealbnd(eqtot_t1fix,[0.1,.1],lb,ub,options);
%wopt = wt2_opt(1)
%t2opt = (wt2_opt(2))


sols = cell(1,Nsweep);

%figure;

parfor n=1:Nsweep
    t1fix = tauvals(n);

    eqtot_t1fix = @(wt2) eqtot((wt2(1)),t1fix,(wt2(2)));


    lb=[0 0];
    ub =[10 10];

    % options.StepTolerance=1e-8;
    %options.ConstraintTolerance=1e-8;

    problem = createOptimProblem('fmincon',...
        'objective',eqtot_t1fix,...
        'x0',[1 1],...
        'lb',lb,'ub',ub);%,'options',options);

    ms = MultiStart;
    nStarts = 500;
    [x,fval,eflag,output,manymins] = run(ms,problem,nStarts);
    exitFlags = [manymins.Exitflag];

    Fvals = [manymins.Fval];

    vals_to_take = find((Fvals<-10)&(exitFlags>0));
    Xvals = [manymins.X];
    ll = length(Xvals)/2;
    XX = reshape([manymins.X],[2,ll]);
    Xaccept = XX(:,vals_to_take);

    %     for nn=1:length(vals_to_take)
    %         tau1plot = t1fix;
    %         tau2plot = Xaccept(2,nn);
    %         wplot = Xaccept(1,nn);
    %         scatter(tau1plot,tau2plot,25,25*wplot,'blue');
    %         hold on;
    %     end
    sols{n} = Xaccept;

end

scatter_all = [];

for n = 1:Nsweep
    Xn = sols{n};
    Xn = unique(round(Xn,3)','rows')';
    sizeXn = size(Xn);
    nsols = sizeXn(2);
    tau1rep = repmat(tauvals(n),[nsols,1]);
    addToMat = [tau1rep';Xn];
    scatter_all = [scatter_all,addToMat];
end


save('boundary_sweepvals');


figure;
dist_thresh=1;
scatter3(scatter_all(1,:),scatter_all(3,:),scatter_all(2,:),25,'filled');
lines=connect_points(scatter_all, dist_thresh);
hold on;
for l = 1:length(lines)
    linel = lines{l};
    pts_on_line = scatter_all(:,linel);
    plot3(pts_on_line(1,:),pts_on_line(3,:),pts_on_line(2,:));
end
