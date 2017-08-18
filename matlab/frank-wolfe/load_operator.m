function op = load_operator(name,param,normalize)

% load_operator - load operator and their derivatives.
%
%   E(x,a) = 1/(2*lambda)*|Phi(x)*a - y|^2 + |a|_1
%   nabla_a E(x,a) = 1/lambda * Phi(x)' * (Phi(x)*a-y) + 1
%                  = 1/lambda * [ C(x,x)*a-C(x,x0)*a0 ] + 1
%   nabla_x E(x,a) = 1/lambda * diag(a)*Phi'(x)' * (Phi(x)*a-y)
%                  = 1/lambda * a .* [ CD(x,x)*a-CD(x,x0)*a0 ]
%
%   Copyright (c) 2016 Gabriel Peyre

% helpers
mynorm = @(x)norm(x(:));
dotp = @(u,v)real( (u(:)') * v(:) );
% for 2D setups
X1 = @(x)x(1:end/2);
X2 = @(x)x(end/2+1:end);

if nargin<3
    normalize = 1;
end

CDD = @(u,x)[]; % not implemented for 2D setting
        
op.z0 = [0;0]; % center point
op.xlim = [-1;1]; op.ylim = [-1;1];

d = 1; % default -- 1D
switch name
    case 'ideal'
        fc = param;
        op.Phi  = @(x)exp( 2i*pi * (-fc:fc)' * x(:)'  );
        op.PhiD = @(x)diag(2i*pi .* (-fc:fc)') * op.Phi(x);
        op.C  = @(u,x)op.Phi(u)' *op.Phi(x);
        op.CD = @(u,x)op.PhiD(u)'*op.Phi(x);
        op.CDD = @(u,x)op.PhiD(u)'*op.PhiD(x);
    case 'ideal-int' % integrated version
        fc = param;
        % pseudo inverse of differentiation
        h = 1 ./ ( 2i*pi .* (-fc:fc)' ); h(fc+1) = 0;
        op.Phi  = @(x)diag(h) * exp( 2i*pi * (-fc:fc)' * x(:)'  );
        op.PhiD = @(x)exp( 2i*pi * (-fc:fc)' * x(:)'  );
        op.C  = @(u,x)op.Phi(u)' *op.Phi(x);
        op.CD = @(u,x)op.PhiD(u)'*op.Phi(x);
        op.CDD = @(u,x)op.PhiD(u)'*op.PhiD(x);
        
    case 'gaussian'
        s = param;
        op.Phi  = @(x)0;
        op.PhiD = @(x)0;
        C  = @(u,x)exp(-(u-x).^2/(2*s^2));
        CD = @(u,x)-(u-x)/(s^2) .* C(u,x);
        CDD = @(u,x)1/s^2 * ( 1 - (u-x).^2/(s^2) ) .* C(u,x);
        
    case 'laplace'
        P = 1024;
        t = 10*(0:P-1)'/P;
        op.Phi  = @(x)exp( -t(:) * x(:)'  );   % should not be used -- discretized
        op.PhiD = @(x)diag(-t) * op.Phi(x);     % should not be used -- discretized     
        C   = @(u,x)1./(u+x); 
        CD  = @(u,x)-1./(u+x).^2;  
        CDD = @(u,x)2./(u+x).^3;             
    case 'laplace-normalized'
        if nargin<2
            param=0;
        end
        a = param; % 0=non-integrable, 0 excluded
        op.Phi  = @(x)0;   % should not be used -- discretized
        op.PhiD = @(x)0;     % should not be used -- discretized
        C   = @(x,y)2*sqrt((x+a).*(y+a)) ./ (x+y+2*a);
        % (y (-x + y))/(Sqrt[x y] (x + y)^2)
        CD  = @(x,y)y.*(y-x)./( sqrt(x.*y) .* (x+y).^2 );
        CD  = @(x,y)CD(x+a,y+a); % take into account the weights
        CDD = @(x,y)( 6*x.*y - x.^2 - y.^2 ) ./ ( 2 * sqrt(x.*y).*( x+y ).^3 );
        CDD = @(x,y)CDD(x+a,y+a); % take into account the weights    
    case 'laplace-int' % integrated version
        op.Phi  = @(x)NaN;   % should not be used -- discretized
        op.PhiD = @(x)NaN;     % should not be used -- discretized     
        C   = @(x,y)x.*(log(x+y)-1) - x.*(log(x)-1) + ...
                    y.*(log(x+y)-1) - y.*(log(y)-1); 
        CD  = @(x,y)log(1+y./x);  
        CDD = @(x,y)1./(x+y);  
    case 'laplace-normalized-int' % integrated version
        op.Phi  = @(x)NaN;   % should not be used -- discretized
        op.PhiD = @(x)NaN;     % should not be used -- discretized     
        C   = @(x,y)-x.^2 .* atan(sqrt(y./x)) + x.*sqrt(x.*y) ...
                    -y.^2 .* atan(sqrt(x./y)) + y.*sqrt(y.*x);
        CD  = @(x,y)2*sqrt(x.*y) - 2.*x.*atan(sqrt(y./x));  
        CDD = @(x,y)sqrt(x.*y)./(x+y);  
        
    case 'gaussian-2d'
        s = param;
        if 0
            %% EXPLICIT FORMULA %%
            C  = @(u1,u2,x1,x2)exp(-( (u1-x1).^2 + (u2-x2).^2 ) /(2*s^2));
            CD1 = @(u1,u2,x1,x2)-(u1-x1)/(s^2) .* C(u1,u2,x1,x2);
            CD2 = @(u1,u2,x1,x2)-(u2-x2)/(s^2) .* C(u1,u2,x1,x2);
        else
            %% USE SYMBOLIC CALCULUS %%
            syms v1 v2 y1 y2;
            CC(v1,v2,y1,y2) = exp(-( (v1-y1)^2 + (v2-y2)^2 ) / (2*s^2));
        end        
        d = 2;

    case 'neuro-like-disc'
        op.z0 = [.4; .3];
        op.xlim = [-1;1];
        op.ylim = xlim;   
        syms v1 v2 y1 y2;     
        C0(v1,v2,y1,y2) = 1-(v1^2+v2^2)*(y1^2+y2^2);
        C1(v1,v2,y1,y2) = (1-y1^2-y2^2)*(1-v1^2-v2^2)*( (1-v1*y1-v2*y2)^2+(v2*y1-v1*y2)^2 );
        CC(v1,v2,y1,y2) = 2*pi*C0(v1,v2,y1,y2)/C1(v1,v2,y1,y2);
        d=2;
        domain = @(x,y)x.^2+y.^2<=.99;
        
    case 'gmixture'
        %% evaluation of x=mean and y=cov of 1D Gaussians         
        op.xlim = [-3;3];
        op.ylim = [.1;4].^2;  
        op.z0 = [0;8];
        syms v1 v2 y1 y2; 
        CC(v1,v2,y1,y2) = 1/sqrt(v2+y2) * exp(- ( (v1-y1)^2 )/(2*( v2+y2 )) );
        d = 2;
        
        
    case 'gmixture2'
        syms v1 v2 y1 y2; 
        CC(v1,v2,y1,y2) = 1/sqrt(v2^2+y2^2) * exp(- ( (v1-y1)^2 )/(2*( v2^2+y2^2 )) );
        op.xlim = [-3;3];
        op.ylim = [.1;4];
        op.z0 = [0;2];
        d = 2;
        
        
    case 'gausslaplace'
        % Gaussian in X and Laplace in Y.
        s = 1; % width of the gaussian
        syms v1 v2 y1 y2; 
        CC(v1,v2,y1,y2) = exp(- (v1-y1)^2 /(2*s^2) ) / (v2+y2);
        d = 2;
        op.z0 = [0;2.5];
        op.xlim = [-3;3];
        op.ylim = [.1 5];
        
    otherwise
        error('Unknown');
end

if exist('CC')
    
    if normalize
        % normalize the columns of Phi in L^2 norm, so covariance is 1 along
        % the diagonal
        CC(v1,v2,y1,y2) = CC(v1,v2,y1,y2)/sqrt( CC(v1,v2,v1,v2) * CC(y1,y2,y1,y2) );
    end

    % USE SYMBOLIC CALCULUS
    % turn formula into evaluable strings
    Cs = fast_formula_2(CC);
    Cs1 = fast_formula_2(diff(CC,v1));
    Cs2 = fast_formula_2(diff(CC,v2));
    %
    C  = @(u1,u2,x1,x2)eval_formula(Cs,u1,u2,x1,x2);
    CD1 = @(u1,u2,x1,x2)eval_formula(Cs1,u1,u2,x1,x2);
    CD2 = @(u1,u2,x1,x2)eval_formula(Cs2,u1,u2,x1,x2);
end

if exist('CDD')
    switch d
        case 1
            remaper = @(f,x,y)f( x(:)*ones(1,length(y)), ones(length(x),1)*y(:)' );
            op.C   = @(x,y)remaper(C,x,y);
            op.CD  = @(x,y)remaper(CD,x,y);
            op.CDD = @(x,y)remaper(CDD,x,y);
        case 2
            op.remaper = @(f,x1,x2,y1,y2)f( x1(:)*ones(1,length(y1)), x2(:)*ones(1,length(y2)), ...
                ones(length(x1),1)*y1(:)', ones(length(x2),1)*y2(:)' );
            op.C   = @(x,y)op.remaper(C,X1(x),X2(x),X1(y),X2(y));
            op.CD  = @(x,y)[op.remaper(CD1,X1(x),X2(x),X1(y),X2(y)); ...
                op.remaper(CD2,X1(x),X2(x),X1(y),X2(y))];
            op.CDD = @(x,y)[];
        otherwise
            error('Not implemented');
    end

op.param = param;

% eta(u)=1/lambda*Phi_u^*( y-Phi_x m_ax )
op.eta  = @(u,x,a,x0,a0,lambda)-1/lambda * real( op.C(u,x)*a -op.C(u,x0)*a0 );
op.etaD = @(u,x,a,x0,a0,lambda)-1/lambda * real( op.CD(u,x)*a-op.CD(u,x0)*a0 );
%%% op.loss = @(x,a,x0,a0)mynorm(op.Phi(x)*a - op.Phi(x0)*a0)^2;
% |Phi(x)*a - op.Phi(x0)*a0|^2 =
% <C(x,x)*a,a> + <C(x0,x0)*a0,a0> - <C(x0,x)*a,a0> - <C(x0,x)*a,a0>
op.loss  = @(x,a,x0,a0)dotp(op.C(x,x)*a,a) ...
    + dotp( op.C(x0,x0)*a0,a0 ) ...
    - dotp( op.C(x0,x)*a,a0 ) ...
    - dotp( op.C(x0,x)*a,a0 );
% TODO: check this is correct in 2D...
op.nablaEx = @(x,a,x0,a0,lambda)1/lambda * repmat(a,[d 1]) .* ( op.CD(x,x)*a-op.CD(x,x0)*a0 );
op.nablaEa = @(x,a,x0,a0,lambda)1/lambda * ( op.C(x,x)*a-op.C(x,x0)*a0 ) + sign(a);
op.E = @(x,a,x0,a0,lambda)1/(2*lambda)*op.loss(x,a,x0,a0) + norm(a,1);
% Vanishing certificate -- only implemented in 1D... TODO
op.etaV = @(u,x0,a0)real( [op.C(u,x0) op.CD(x0,u)'] * ...
        ( [op.C(x0,x0),  op.CD(x0,x0)'; ...
         op.CD(x0,x0) op.CDD(x0,x0)] \ ...
         [sign(a0); a0*0] ) );
     
% energy of a single dirac
op.E_single = @(x,x0,lambda)-(op.C(x,x0)-lambda).^2 ./ (2*diag(op.C(x,x))) + op.C(x0,x0)/2;

% optimal height, assuming same sign
op.Optimal_a = @(x,x0,a0,lambda)op.C(x0,x0)\(op.C(x0,x0)*a0-lambda*sign(a0));

end