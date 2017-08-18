function EtaW = compute_etaw(C,P,z0)

% compute_etaw - compute etaW certificate in 1-D and 2-D
%
%   EtaW = compute_etaw(C,P,z0);
%
%   C(x,y,x1,y1) should be a symbolic function.
%   EtaW(x,y) is a symbolic function.
%   Each P(i) should be a symbolic polynomial in (x,y).
%   z0 is the center point (default is [0;0])
%
%   For 1-D problems, simply make the above functions not depend on (y,y1).
%
%   Copyright (c) 2017 Clarice Poon and Gabriel Peyre

syms x y x1 y1;

if nargin<3
    z0 = [0 0]';
end

R = length(P);

% build the Gamma matrix
Gamma = zeros(R,R);
s = 0;  G = {};
fprintf('Building Gamma: ');
for a=1:R
    % differentiate on (x,y)
    G{a} = poly_deriv(C,P(a));
    for b=a:R
        s = s+1; progressbar(s,R*(R+1)/2);
        % differentiate on (x1,y1), so we switch the variable
        H = poly_deriv(G{a}(x1,y1,x,y),P(b));
        % evaluate at (x,y,x1,y1)=(z0,z0)
        Gamma(a,b) = double( subs(H,{x,y,x1,y1},{z0(1),z0(2),z0(1),z0(2)}) ); 
        Gamma(b,a) = Gamma(a,b); % Gamma is a symmetric matrix.
    end
end

if not(issymmetric(Gamma))
    error('Gamma is not symmetric');
end
if cond(Gamma)>1e6
    warning('Gamma is singular or close to singular.');
end

% coefficients
U = pinv(Gamma)*[1;zeros(R-1,1)];
% build etaW
etaW = 0; 
fprintf('Building etaW: ');
for a=1:R
    progressbar(a,R);
    etaW = etaW + U(a)*subs(G{a},{x y x1 y1},{z0(1) z0(2) x y});
end
EtaW = symfun(etaW,[x y]);

end