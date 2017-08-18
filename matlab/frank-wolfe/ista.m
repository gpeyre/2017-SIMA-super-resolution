function [a,Err,A] = ista(op, x,a,x0,a0, lambda,options)

% ista - iterative soft thresholding
%
%   [a,Err] = ista(op, x,a,x0,a0, lambda,options);
%  
%   min_a 1/2*|Phi_x*a-Phi_x0*a0|^2 + lambda*|a|
%
%   options.niter is #iter.
%   x is the grid on which the solution is discretized
%   a is the initialization (can be 0)
%
%   Copyright (c) 2016 Gabriel Peyre

pos_constr = getoptions(options, 'pos_constr', 0);

if pos_constr==1
    % impose positivity
    ProxL1 = @(a,gamma)max(a-gamma,0);
else
    % plain L1
    ProxL1 = @(a,gamma)max(abs(a)-gamma,0) .* sign(a);
end

options.null = 0;
niter = getoptions(options, 'niter', 100);
tau = getoptions(options, 'tau', -1);
verb = getoptions(options, 'verb', 0);
if tau<0
    tau = 1.5/norm(op.C(x,x));
end

Cxx  = op.C(x,x);
Cx0x = op.C(x0,x);
Cx0x0 = op.C(x0,x0);
r = op.C(x,x0)*a0;
dotp = @(u,v)real( (u(:)') * v(:) );
loss = @(a)dotp(Cxx*a,a)+dotp(Cx0x0*a0,a0)-2*dotp(Cx0x*a,a0);

Err = []; A = [];
for j=1:niter
    if verb==1
        progressbar(j,niter);
    end
    G = real( Cxx*a-r );
    a = ProxL1( a-tau*G, lambda*tau );
    if nargout>2
        A(:,end+1) = a;
    end
    % Err(end+1) = op.E(x,a,x0,a0,lambda);
    Err(end+1) = 1/(2*lambda)*loss(a)+norm(a,1);
end

end