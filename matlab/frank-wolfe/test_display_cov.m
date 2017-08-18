%%
% Display of covariance functions.

test = 'convolution';
test = 'laplace';

switch test
    case 'laplace'
        name = {'laplace-normalized', ...
            'laplace', ...
            'laplace-int', ...
            'laplace-normalized-int', ...
            };
        dom = [.1 5];
        param_op = 0;
    case 'convolution'
        name = {'ideal-int', ...
            'ideal', ...
            };
        dom = [0 1];
        param_op = 8; % =fc
end

rep = 'results/covariances-display/';
if not(exist(rep))
    mkdir(rep);
end

x = linspace(dom(1),dom(2),512);

for k=1:length(name)
    op = load_operator(name{k},param_op);
    
    %% display covariance 
    clf; hold on;
    c = real( op.C(x,x) );
    imagesc(x, x, c );
    contour(x,x,c, 20, 'k');
    colormap parula(256);
    box on; axis image; 
    drawnow;
    saveas(gcf, [rep name{k} '-cov.png'], 'png');
    
    %% display etaV
    X0 = {[.5], [.3 .7]', [.25 .5 .7]', [.25 .5 .7]', [.25 .5 .7]'}; 
    A0 = {[1], [1 1]', [1 1 1]', [1 1 -1]', [1 -1 1]'};
    for r=1:length(X0)
        x0 = X0{r}; a0 = A0{r};
        x0 = dom(1) + x0*(dom(2)-dom(1));
        etaV = op.etaV(x, x0, a0);
        clf; hold on;
        plot(x, etaV);
        stem(x0, a0, 'r:.', 'MarkerSize', 10);
        plot([min(x), max(x)], [1 1], 'k--');
        plot([min(x), max(x)],-[1 1], 'k--');
        axis([min(x), max(x), -1.05, 1.05]);
        box on;
        drawnow;
        saveas(gcf, [rep name{k} '-etaV-' num2str(r) '.png'], 'png');
        if length(x0)>1 %% ZOOM
            vmin = min(etaV(x>x0(1) & x<x0(end)));
            axis([x0(1) x0(end) vmin 1]);
            saveas(gcf, [rep name{k} '-etaV-' num2str(r) '-zoom.png'], 'png');
        end
    end

end