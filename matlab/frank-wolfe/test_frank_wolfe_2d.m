%%
% Test for Frank-Wolfe on BLASSO in 2D

addpath('toolbox/');
addpath('minConf/');

name = 'gaussian-2d';
name = 'gausslaplace';
name = 'gmixture2';
name = 'neuro-like-disc';
name = 'gmixture';
param_op = .3; % width of the Gaussian

rep = ['results/fw-2d/' name '/'];
[~,~] = mkdir(rep);

%%
% Helpers.

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
myplot = @(x,y,c)area(x, y, 'FaceColor', c, 'EdgeColor', c, 'FaceAlpha', 1, 'EdgeAlpha', 1);
lw = 2;
X1 = @(x)x(1:end/2);
X2 = @(x)x(end/2+1:end);

%%
% Parameters.

% load the covariance operators
op = load_operator(name,param_op);

% number of spikes
N = 2;
N = 3;
% gris seach space for spikes addition
P = 256;
xu = linspace(op.xlim(1),op.xlim(2),P);
yu = linspace(op.ylim(1),op.ylim(2),P);
[Yu,Xu] = meshgrid(yu,xu);
U = [Xu(:);Yu(:)];

% input measure
test_type = 'randn';
test_type = 'equi';
switch test_type
    case 'randn'
        x0c = randn(N,2);
        x0c = x0c - repmat(mean(x0c,1), [N 1]);
        x0c = x0c / max( sqrt(sum(x0c.^2,2)) );
    case 'equi'
        t = pi/8 + 2*pi*(0:N-1)'/N;
        x0c = [cos(t),sin(t)];       
end

x0c(:,1) = x0c(:,1)*diff(op.xlim)/2;
x0c(:,2) = x0c(:,2)*diff(op.ylim)/2;

% scale positions
delta = .4; % spacing
x0 = repmat( op.z0(:)', [N 1] ) + delta*x0c;
x0 = x0(:);
% amplitude
a0 = ones(N,1);

%%
% Observations.

y = real( op.C(U,x0)*a0 ); y = y/max(y);
y = reshape(y,[P P]);
clf; hold on;
imagesc(xu,yu,y'); 
colormap parula(256);
axis tight; axis square; 
plot(X1(x0), X2(x0), 'r.', 'MarkerSize', 20);

%%
% Regularization parameter

lambda = .01*max(y(:));
lambda = .0001*max(y(:));

%% 
% Resolution using F-W.

options.pos_constr = 1; % impose positivity or not
options.fw_grid = U;        % grid to intialize the seed search
options.niter_ista = 300;   % sup-iterations when using ista
options.niter_bfgs = 600;   % sup-iterations when using bfgs
options.niter = 20*N;        %  max #iterations of FW
options.gridding_refine = 1;    % refine the grid search for the FW step by BFGS
options.tol = 1e-8;         % tolerance on l^inf norm of certificate for stopping
options.fw_update = 'non-convex';   % way to update (position,amplitude)
options.display = 1;
options.bfgs_solver = 'minconf';
options.bfgs_solver = 'hanso';
[x,a,X,A] = frank_wolfe(op, x0, a0,lambda,options);

%%
% Display computed solution.


% eta_lambda, should approximation eta0
eta1 = op.eta(U,x,a,x0,a0,lambda);
eta1 = reshape(eta1, [P P]);

% display
clf;
hold on;
imagesc(xu,yu,eta1'); 
contour(xu,yu,eta1', 20, 'k');
colormap parula(256);
plot_measure_2d(x,a,'b',50);
plot(X1(x0), X2(x0), 'ro', 'MarkerSize', 5);
set(gca, 'FontSize', 15);
axis tight; axis square; box on;


%%
% Evolution with lambda.

%%% WARNING: need to tune a bit this for each plot %%%
options.niter = 4*N; % impose strict sparsity, but might not reach convergence 
options.display = 0;
K = 40;
lambda_max = max(y(:))*.2;
lambda_list = linspace(0.005,1,K)*lambda_max;
X = {}; A = {};
for k=1:K
    progressbar(k,K);
    [X{end+1},A{end+1}] = frank_wolfe(op, x0, a0,lambda_list(k),options);
end

% remove those with too little #spikes
I = find(cellfun(@length,X)>=2*N);
X = {X{I}}; A = {A{I}};
K = length(I);


% display
clf;
hold on;
imagesc(xu,yu,y'); 
contour(xu,yu,y', 20, 'k');
colormap parula(256);
for k=1:K
    t = (k-1)/(K-1);
    xk = reshape(X{k},[length(X{k})/2 2]);
    plot(xk(:,1), xk(:,2), '.', 'color', [t,0,1-t], 'MarkerSize', 20);
end
set(gca, 'FontSize', 15);
axis tight; axis square; box on;
plot(X1(x0), X2(x0), 'k.', 'MarkerSize', 30);
saveas(gcf, [rep name 'delta' num2str(delta) '.png'], 'png');
zoom(.9/delta); box on;
saveas(gcf, [rep name 'delta' num2str(delta) '-zoom.png'], 'png');