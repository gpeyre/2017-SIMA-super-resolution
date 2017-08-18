function [x,a,R] = noncvx_sparse_spikes(op,lambda, x0,a0, x1,a1, options)

% [xf,af] = noncvx_sparse_spikes(op,x0,a0, x1,a1)
%
%   Solve 
%       min_{x,a} 1/(2*lambda)|Phi_x0*a0-Phi_x*a|^2 + |a|_1
%   using L-BFGS with initialization (x1,a1).
%
%   Set options.niter for BFGS #iterations.
%
%   Copyright (c) 2017 Gabriel Peyre

options.null = 0;
bfgs_solver = getoptions(options, 'bfgs_solver', 'hanso');

X = @(z)z(1:length(x1));
A = @(z)z(length(x1)+1:end);
XA = @(x,a)[x(:);a(:)];

z1 = XA(x1,a1);
Ebfgs = @(z)op.E(X(z),A(z),x0,a0,lambda);
nablaE = @(z)real( XA( op.nablaEx(X(z),A(z),x0,a0,lambda), op.nablaEa(X(z),A(z),x0,a0,lambda) ) );
% callback for L-BFGS
nablaEbfgs = @(z)deal(Ebfgs(z),nablaE(z));

switch bfgs_solver
    case 'hanso'
        options.report = @(z,v)v;
        [z, R, info] = perform_bfgs(nablaEbfgs, z1, options);
    case 'minconf'
        N = length(a1);
        LB = [ones(N,1)*op.xlim(1);ones(N,1)*op.ylim(1);zeros(N,1)-Inf];
        UB = [ones(N,1)*op.xlim(2);ones(N,1)*op.ylim(2);zeros(N,1)+Inf];
        opt.verbose = 0;
        opt.maxIter = getoptions(options, 'niter_bfgs', 100);
        opt.optTol = 1e-12;
        [z,R] = minConf_TMP(nablaEbfgs,z1,LB,UB,opt);        
end
x = X(z); a = A(z);

end