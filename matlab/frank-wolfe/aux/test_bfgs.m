%%
% test for BFGS for non convex sparsity.

name = 'laplace';
name = 'ideal';
name = 'laplace-normalized';

addpath('toolbox/');

% helpers
mynorm = @(x)norm(x(:));


% just display
P = 512*4;
u = (0:P-1)'/P;

fc = 10; 
N = 2; % # diracs


% input diracs
switch name
    case 'ideal'
        param_op = fc;
        delta = .3;
        x0 = [.5-delta/2; .5+delta/2];
        a0 = [1 1]';
    case {'laplace' 'laplace-normalized'}
        param_op = 0.1; % prevent from problem with spike at 0.
        switch N
            case 2
                x0 = [1 5]'; a0 = [.9;1.1];
                x1 = [.9;5.1]; a1 = [1;1]; % init
                u = rescale(u,0,7);
            case 1
                x0 = [1]; a0 = [1];
                x1 = x0; a1 = a0; 
                u = rescale(u,0,3);
        end
end

% only used for 'ideal'
op = load_operator(name,param_op);

%%
% Minimal norm certificate.

etaW = op.etaV(u,x0,a0);
clf; hold on;
plot(u, etaW); 
stem(x0, a0, 'r:.', 'MarkerSize', 10);
plot([min(u) max(u)], [1 1], 'k--');
axis([min(u), max(u), -.1, 1.1]);
box on;

%% 
% Non-convex resolution 

lambda = .1; 

% option for BFGS
options.niter = 100;
[xf,af,R] = noncvx_sparse_spikes(op,lambda, x0,a0, x1,a1, options);

clf; 
plot(R); 
axis tight;

% compare against energy at seed point
a2 = op.Optimal_a(x0,x0,a0,lambda);
fprintf('Initial/Optimized-1=%.3f (should be>=0)\n', op.E(x0,a2,x0,a0,lambda) / op.E(xf,af,x0,a0,lambda) - 1);

% should be <1 and interpolating to certify global minimizer 
eta1 = op.eta(u,x1,a1,x0,a0,lambda);
etaf = op.eta(u,xf,af,x0,a0,lambda);

clf; hold on;
plot(u, eta1, 'b:');
plot(u, etaf, 'r');
stem(x1, a1, 'b:.', 'MarkerSize', 10);
stem(xf, af, 'r.', 'MarkerSize', 10);
plot([min(u) max(u)], [1 1], 'k--');
axis([min(u) max(u) -.4 1.4]);
box on;

%%
% display energy

plot(u, op.E_single(u,x0,lambda) );

t = linspace(.05,1,100);
e = arrayfun(@(t)op.E(t*x0,a0,x0,a0,lambda), t);
clf;
plot(t,e);



