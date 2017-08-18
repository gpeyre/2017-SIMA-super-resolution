tic
addpath('toolbox/');

name = 'neuro-like-disc';

rep = ['results/examples/' name '/'];
[~,~] = mkdir(rep);

% number of points
N = 1;
% N = 6;
N = 2;
% N = 3;

align = 1;

randn('state', 321); % comment this to avoid generating same points

%points
if align    
    theta = pi*rand();
    Z = 10*(1:N)'*[cos(theta) sin(theta)];
else
    Z = round(randn(N,2)*100);
end
%center
z0 = [0.4,0.3];
d=2;

%%
%%
% Test for convergence of eta_V
xlim = [-1,1];
ylim = [-1,1];
tMax = .25;
Nt = 5;
for k=1:Nt
    %     progressbar(k,Nt);
    t = tMax*k/Nt;
    % compute interpolation points.
    if N~=1
    Zc = Z-repmat(mean(Z,1), [size(Z,1) 1]);
    Zc = Zc/max(Zc(:));
    x0 = repmat(z0(:)',[N 1]) + t*Zc;
    else
        x0 = z0;
    end
    % computes eta_V
    EtaV = compute_etav_neuro_disc(x0);
    
    % display
    clf;
    display_eta(EtaV,d,Z,z0,x0,xlim,ylim);
    drawnow;
    saveas(gcf, [rep name '-etav-N' num2str(N) '-align' num2str(align) '-' num2str(k) '.png'], 'png');

    if N==1
        break;
    end
end
toc

%%
