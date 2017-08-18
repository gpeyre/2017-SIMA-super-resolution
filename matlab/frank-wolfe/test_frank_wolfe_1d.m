%%
% test for Frank-Wolfe on BLASSO

% E(m) = 1/(2*lambda)*|Phi*m - y|^2 + |m|_1
%   partial E(m) = 1/(2*lambda)*Phi^*( Phi*m - y ) + |m|_1


name = 'laplace';
name = 'laplace-normalized';
name = 'ideal';
name = 'gaussian';

rep = 'results/fw/';
[~,~] = mkdir(rep);

N=8;
N=4;
N=1;
N=3;

addpath('toolbox/');

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
myplot = @(x,y,c)area(x, y, 'FaceColor', c, 'EdgeColor', c, 'FaceAlpha', 1, 'EdgeAlpha', 1);
lw = 2;

% gris seach space for spikes addition
P = 512*32;
u = (0:P-1)'/P;

% input diracs
switch name
    case 'old'
        fc = 10;
        param_op = fc;
        delta = .3;
        x0 = [.5-delta/2; .5+delta/2];
        a0 = [1 1]';
    case {'ideal' 'gaussian'}
        if strcmp(name, 'gaussian')
            param_op = .1;
        else
            param_op = 10; % cutoff
        end
        switch N
            case 1
                x0 = .5
                a0 = 1;
            case 2
                x0 = [.5-delta/2; .5+delta/2];
                a0 = [1 1]';
            case 3
                x0 = [.3 .5 .8];
                a0 = [.8 1 -.8]';
            case 4
                x0 = [.15 .35 .6 .9];
                a0 = [-.8 .8 1 -.8]';
            otherwise
                x0 = sort(rand(N,1)); 
                a0 = sign(randn(N,1)) .* (.5+.5*rand(N,1));
        end
    case {'laplace' 'laplace-normalized'}
        param_op = 0.1; % prevent from problem with spike at 0.
        switch N
            case 1
                x0 = [1]; a0 = [1];
                u = rescale(u,0,3);
            case 2
                x0 = [1 5]'; a0 = [.9;1.1];
                u = rescale(u,0,7);
            case 3
                x0 = [1 3 5]'; a0 = [1.2;.8;1.1];
                u = rescale(u,0,7);
        end
end

op = load_operator(name,param_op);

%%
% Observations.

y = real( op.C(u,x0)*a0 ); y = y/max(y);
clf; hold on;
plot(u,y, 'LineWidth', lw);
stem(x0, a0, 'r:.', 'MarkerSize', 25, 'LineWidth', lw);
axis tight;
set(gca, 'XTick', [], 'YTick', []); box on;
SetAR(1/3);    
box on;
saveas(gcf, [rep name '-N' num2str(N) '-observ.eps'], 'epsc');


%%
% Minimal norm certificate.

etaV =op.etaV(u,x0,a0);
clf; hold on;
plot(u, etaV, 'LineWidth', lw); 
stem(x0, a0, 'r:.', 'MarkerSize', 25, 'LineWidth', lw);
plot([min(u) max(u)], [1 1], 'k--');
plot([min(u) max(u)],-[1 1], 'k--');
axis([min(u) max(u) min(etaV)-.05 max(etaV)+.05]); 
set(gca, 'XTick', [], 'YTick', []); box on;
SetAR(1/3);    
box on;
saveas(gcf, [rep name '-N' num2str(N) '-etaV.eps'], 'epsc');

%%
% Regularization parameter (small on purpose!)

lambda = .01;

%% 
% Resolution using ISTA on a discretization grid.

do_ista = 0;

if do_ista
    
Q = 256*2;
z = (0:Q-1)'/Q; 

T = 12000*20; % #iter
opt.niter = T;
opt.tau = 1.2/norm(op.C(z,z)); % slow convergence
opt.verb = 1;
[a1,Err,A1] = ista(op, z,z*0,x0,a0, lambda,opt);

% re-sample position to obtain smoother video
u = 2;
R = sqrt( sum( abs(A1(:,1:end-1)-A1(:,2:end)).^u ) );
R = [0 cumsum(R)]; R = R/R(end);
K = 60; % #frame in video
m = interp1(R(:)',1:T, linspace(0,1,K)); 
% plot and generate video
F = [];
s = 0;
for k=1:K
    t = (k-1)/(K-1);
    % interp to get value
    u = floor(m(k)); v = ceil(m(k)); s = m(k)-u;
    a1 = (1-s)*A1(:,u) + s*A1(:,v);    
    %
    clf; hold on;
    c = [t 0 1-t];
    myplot(z,a1,c);
%    h = bar(z,a1,'FaceColor', c, 'EdgeColor', c); axis tight;
    % set(gca, 'FontSize', 20);
    axis([0 1 min(A1(:)) max(A1(:))]);
    set(gca, 'XTick', [], 'YTick', []); box on;
    SetAR(1/3);
    drawnow;
    f = getframe;
    F(:,:,:,k) = f.cdata;
end
% find colormap
A = permute(F, [1 2 4 3]);
A = reshape(A, [size(A,1) size(A,2)*size(A,3) 3]);
[~,map] = rgb2ind(uint8(A),254,'nodither');
map(end+1,:) = 0; map(end+1,:) = 1;
% convert
im = [];
for s=1:size(F,4);
    im(:,:,1,s) = rgb2ind(uint8(F(:,:,:,s)),map,'nodither');
end
% save
imwrite(im+1,map,[rep name '-N' num2str(N) '-ista.gif'], ...
        'DelayTime',0,'LoopCount',inf);
    
end

%% 
% Resolution using F-W.

options.fw_grid = u;        % grid to intialize the seed search
options.niter_ista = 300;   % sup-iterations when using ista
options.niter_bfgs = 100;   % sup-iterations when using bfgs
options.niter = 100;        % iterations of FW
% way to update (position,amplitude)
options.fw_update = 'ista';
options.fw_update = 'non-convex';

[x,a,X,A] = frank_wolfe(op, x0, a0,lambda,options);


nor = @(x)x/max(abs(x));
X1 = {.5, X{:}};
A1 = {0, A{:}};
ms = 25; lw = 2; 
for k=1:length(A1)
    t = (k-1)/(length(A1)-1);
    x = X1{k}; a = A1{k};
    eta1 = op.eta(u,x,a,x0,a0,lambda);
    V = max(abs(eta1));
    clf; hold on;
    plot(u, eta1, 'color', [t 0 1-t], 'LineWidth', lw);
    stem(x0, nor(a0)*V, 'k:.', 'MarkerSize', ms, 'LineWidth', lw);
    if max(abs(a(:)))>0
        h = stem(x, nor(a)*V, '.', 'MarkerSize', ms, 'LineWidth', lw);
        set(h, 'color', [t 0 1-t]);
    end
    plot([min(u) max(u)], [1 1], 'k--');
    plot([min(u) max(u)],-[1 1], 'k--');
    axis([min(u), max(u), -V*1.03, V*1.03]);
    set(gca, 'XTick', [], 'YTick', []); box on;
    SetAR(1/3);
    saveas(gcf, [rep name '-N' num2str(N) '-fw-' num2str(k) '.eps'], 'epsc');
    % set(gca, 'FontSize', 20);
end


return;


T = 1000; % #iter
opt.niter = T;
opt.tau = 1.2/norm(op.C(z,z)); % slow convergence
opt.verb = 1;
[a1,Err,A1] = ista(op, z,z*0,x0,a0, lambda,opt);

opt.niter = 2;
opt.tau = 1.2/norm(op.C(z,z)); % slow convergence
[a1,Err,A1] = ista(op, z,a1,x0,a0, lambda,opt);

