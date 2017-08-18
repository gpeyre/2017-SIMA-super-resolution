
addpath('toolbox/');


rep = 'results/gmixtures/';
[~,~] = mkdir(rep);

if not(exist('name'))
name = 'gmixture';
name = 'gmixture2';
name = 'gmixture3';
end

N = 2;
normalize = 1;
[C,d,z0,xlim,ylim] = load_correlation(name,normalize);

%%
% Compute vanishing derivatives constraints

% degree of expansion
K = poly_order(N)+1;

M = 10;
for k=1:M
    % direction
    theta = pi/2*(k-1)/(M-1);
    Z = [z0(1) z0(2);
        z0(1) + cos(theta), z0(2) + sin(theta)];
    Z = round(100*Z);
    % compute
    [M_fac,C_sep, Nodes,max_deg] = TaylorMtx(Z,K);    
    P = C_sep{1};
    EtaW = compute_etaw(C,P,z0);
    % display
    clf;
    display_etaw(EtaW,d,Z,z0,[];xlim,ylim);
    saveas(gcf, [rep name '-N' num2str(N) '-nor' num2str(normalize) '-' num2str(k) '.png'], 'png');
end


