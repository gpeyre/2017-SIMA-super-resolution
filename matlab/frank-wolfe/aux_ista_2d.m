
%% 
% Resolution using ISTA on a discretization grid.

do_ista = 0;
if do_ista
    % discretization grid
    Q = 64;
    z = (0:Q-1)'/Q;
    [Yz,Xz] = meshgrid(z,z);
    Z = [Xz(:);Yz(:)];
    % ISTA solver
    T = 1200; % #iter
    opt.niter = T;
    opt.tau = 1.2/norm(op.C(Z,Z)); % slow convergence
    opt.verb = 1;
    opt.pos_constr = 1; % impose positivity or not
    [a1,Err,A1] = ista(op, Z,zeros(length(Z)/2,1),x0,a0, lambda,opt);
    A1 = reshape(A1, [Q Q size(A1,2)]);
    a1 = reshape(a1,[Q Q]);    
    % video!
    k = 0;
    while true
        k = mod(k,T)+1;
        imagesc(A1(:,:,k)); axis image; axis off;
        colormap parula(256); drawnow;
    end
end