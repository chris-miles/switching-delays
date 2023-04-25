function [y,t,jump_times,states] = run_1stochsim(tau1,tau2,alpha,beta,epsilon,I, K, d, nn, M,Tmax)

states0 = (rand<(alpha/(alpha+beta)))+1; % flip a coin 1 or 2

Q = [-alpha alpha; beta -beta]/epsilon; % 2 state Markov chain for now
[jump_times,states] = simCTMC(Q,Tmax,1,states0);
while length(states)<2
    [jump_times,states] = simCTMC(Q,Tmax,1,states0);
end 
%limP = limitdist(Q);


%% DDE stuff
delays = [tau1  tau2]; % switch between these two delays

maxdelay = max(delays);



dydt = @(t,y,z)  I-d*y(1) - K*(z(1,1)^nn)/(M^nn+z(1,1)^nn); 

p0 = 0; % initial value
history0 = @(t) p0*ones(1,1); % pick (constant) history


tspan0 = [0,jump_times(2)];
sol = dde23(dydt, delays(states(1)), history0, tspan0);

tagg = [-maxdelay,sol.x]'; 
yagg = [history0(-maxdelay),sol.y]';  % add t=-maxdelay as a data point
                            % need to modify if time-dependent initial
                            % history

njumps = length(jump_times);

for n = 2:njumps
    tstart_n = jump_times(n);
    if n==njumps
        tend_n=Tmax;
    else
        tend_n = jump_times(n+1);
    end
    state_n = states(n);
    delay_n = delays(state_n);
    
    interped_history = @(t) interp1(tagg, yagg, t);
    trange_n = [tstart_n, tend_n];

    sol_n = dde23(dydt, delay_n, interped_history, trange_n);
    soln_t=sol_n.x;
    soln_y = sol_n.y;
    tagg=[tagg;soln_t(2:end)'];
    yagg=[yagg;soln_y(:,2:end)'];
end


t = tagg;
y=yagg;
end 