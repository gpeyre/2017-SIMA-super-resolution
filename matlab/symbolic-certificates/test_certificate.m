
addpath('toolbox/');

name = 'isolaplace';
name = 'gaussian1d';
name = 'laplace1d';
name = 'laplace1d-normalized';
name = 'neuro-like';
name = 'gaussian2d';
name = 'neuro-like-1d-3';
name = 'neuro-like-2d';
name = 'covindep';
name = 'neuro-like-box';
name = 'gaussian2d';
name = 'neuro-like-box';
name = 'gmixture';
name = 'gausslaplace';
name = 'neuro-like-disc';


rep = ['results/examples/' name '/'];
[~,~] = mkdir(rep);

% number of points
N = 7; 
N = 3;
N = 1;
N = 2;

% set to 1 to impose that the diagonal of the covariance is 1
normalize = 1;
[C,d,z0,xlim,ylim,opt.domain] = load_correlation(name,normalize);

%%
% Compute vanishing derivatives constraints

switch d
    case 1
        K = 2*N-1;
        % list of polynomials describing the derivative constraints
        clear P;
        for k=0:K
            P(k+1) = x^k;
        end
    case 2
        % degree of expansion
        K = poly_order(N)+1;
        % points
        randn('state', 321); % comment this to avoid generating same points
        Z = round(randn(N,2)*10);
        % Z = [[0 0]; [1 0]]; 
        [M_fac,C_sep, Nodes,max_deg] = TaylorMtx(Z,K);
        % list of polynomials describing the derivative constraints
        P = C_sep{1};
end

EtaW = compute_etaw(C,P,z0);

%%
% Display.

clf; 
display_eta(EtaW,d,Z,z0,[],xlim,ylim);
saveas(gcf, [rep name '-etaw-N' num2str(N) '-nor' num2str(normalize) '.png'], 'png');

%%
% Test for convergence of eta_V

tMax = .5;
Nt = 5;
for k=1:Nt
    progressbar(k,Nt);
    t = tMax*k/Nt;
    % compute interpolation points.
    Zc = Z-repmat(mean(Z,1), [size(Z,1) 1]);
    Zc = Zc/(1e-10 + max(Zc(:)));
    x0 = repmat(z0(:)',[N 1]) + t*Zc;
    % computes eta_V
    EtaV = compute_etav(C,x0);
    % display
    clf;
    display_eta(EtaV,d,Z,z0,x0,xlim,ylim, opt);
    drawnow;
    saveas(gcf, [rep name '-etav-N' num2str(N) '-nor' num2str(normalize) '-' num2str(k) '.png'], 'png');
end