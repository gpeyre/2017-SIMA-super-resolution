function EtaV = compute_etav_neuro_disc(x0)

% compute_etaw - compute eta_V certificate for neuro-like-disc kernel
% varphi(x,y): t\in [0,2pi] -> (|x-cos(t)|^2 + |y-sin(t)|^2)^{-1}
%   EtaW = compute_etav(C,x0);
%
%   EtaV(x,y) is a symbolic function.
%   Each x0(i,:) is a point.
%
%
%   Copyright (c) 2017 Clarice Poon and Gabriel Peyre

syms x y x1 y1;

syms z p p1;

N = size(x0,1);

C = -1i*z/( (p-z)*(p'*z-1)*(p1-z)*(p1'*z-1)  );
Cx = 2*1i*( real(p)*z^2- z^3/2-z/2 )/( (p-z)^2*(p'*z-1)^2*(p1-z)*(p1'*z-1) );
Cx1 = 2*1i*( real(p1)*z^2- z^3/2-z/2 )/( (p-z)*(p'*z-1)*(p1-z)^2*(p1'*z-1)^2 );


Cy = 2*1i*( imag(p)*z^2- z^3/(2*1i)+z/(2*1i) )/( (p-z)^2*(p'*z-1)^2*(p1-z)*(p1'*z-1) );
Cy1 = 2*1i*( imag(p1)*z^2- z^3/(2*1i)+z/(2*1i) )/( (p-z)*(p'*z-1)*(p1-z)^2*(p1'*z-1)^2 );



Cyy1 = -4*1i*z*( imag(p)*z-z^2/(2*1i)+ 1/(2*1i))*(imag(p1)*z-z^2/(2*1i) + 1/(2*1i) )/...
    ((p-z)^2*(z*p'-1)^2*(p1-z)^2*(z*p1'-1)^2);
Cxx1 = -4*1i*z*( real(p)*z-z^2/2 -1/2)*(real(p1)*z-z^2/2 -1/2 )/...
    ((p-z)^2*(z*p'-1)^2*(p1-z)^2*(z*p1'-1)^2);
Cyx1 = -4*1i*z*( imag(p)*z-z^2/(2*1i)+ 1/(2*1i))*(real(p1)*z-z^2/2 -1/2 )/...
    ((p-z)^2*(z*p'-1)^2*(p1-z)^2*(z*p1'-1)^2);



cres_1 = @(F)( cauchy_residue(subs(F,p1,p), p ) );
cres_2 = @(F)( cauchy_residue(F, [p p1] ) );

evalij_1 = @(F,pt) double(subs(cres_1(F), p, pt));
evalij_2 = @(F,pts) double(subs(cres_2(F), [p p1], pts));

normalization = 1/sqrt(subs(cres_1(C), p, x+1i*y) * subs(cres_1(C), p, x1+1i*y1));
dx_norm = diff(normalization,x,1);
dx1_norm = diff(normalization,x1,1);
dxx1_norm = diff(diff(normalization,x,1),x1,1);
dy_norm = diff(normalization,y,1);
dy1_norm = diff(normalization,y1,1);
dyy1_norm = diff(diff(normalization,y,1),y1,1);
dyx1_norm = diff(diff(normalization,y,1),x1,1);


evalij_norm = @(fn, i,j)  subs(fn, {x y x1 y1}, {x0(i,1) x0(i,2) x0(j,1) x0(j,2)});

% build the Gamma matrix
s=0;
Gamma = (zeros(3*N,3*N));
for i=1:N
    for j=1:N
        s = s+1; progressbar(s,N*(N+1)/2);
        pa = x0(i,1)+1i*x0(i,2);
        pb = x0(j,1)+1i*x0(j,2);
        if abs(pa-pb)<1e-10
            pts = pa;
            evalij = evalij_1;
        else
            pts = [pa,pb];
            evalij = evalij_2;
        end
        
        Gamma(i,j) = evalij(C,pts)*evalij_norm(normalization,i,j);
        Gamma(i+N,j+N) = evalij(Cxx1,pts)*evalij_norm(normalization,i,j) + evalij(Cx,pts)*evalij_norm(dx1_norm,i,j) + ...
        + evalij(Cx1,pts)*evalij_norm(dx_norm,i,j) + evalij(C,pts)*evalij_norm(dxx1_norm,i,j); 
        
        Gamma(i+2*N,j+2*N) =  evalij(Cyy1,pts)*evalij_norm(normalization,i,j) + evalij(Cy,pts)*evalij_norm(dy1_norm,i,j)...
        + evalij(Cy1,pts)*evalij_norm(dy_norm,i,j)+ evalij(C,pts)*evalij_norm(dyy1_norm,i,j); 
              
        Gamma(i+N,j) = evalij(Cx,pts)*evalij_norm(normalization,i,j) + evalij(C,pts)*evalij_norm(dx_norm,i,j); 
        Gamma(j,i+N) = Gamma(i+N,j);
        Gamma(i+2*N,j) = evalij(Cy,pts)*evalij_norm(normalization,i,j) + evalij(C,pts)*evalij_norm(dy_norm,i,j); 
        Gamma(j,i+2*N) = Gamma(i+2*N,j);
        Gamma(i+2*N,j+N) = evalij(Cyx1,pts)*evalij_norm(normalization,i,j)+evalij(Cy,pts)*evalij_norm(dx1_norm,i,j)+...
         evalij(Cx1,pts)*evalij_norm(dy_norm,i,j)+evalij(C,pts)*evalij_norm(dyx1_norm,i,j); 
                
        Gamma(j+N,i+2*N) = Gamma(i+2*N,j+N);   
    end
end

if not(issymmetric(Gamma))
    error('Gamma is not symmetric');
end
cnum = cond(Gamma);
if cnum>1e6
    warning(['Gamma is singular or close to singular: ', num2str(cnum)]);
end


% coefficients
U = pinv(Gamma)*[ones(N,1);zeros(2*N,1)];
% build EtaV
subsij0 = @(F,i)subs(F,{x y x1 y1},{x0(i,1) x0(i,2) x y});
subsij = @(F,i) subs(cres_2(F), [p p1], [x0(i,1)+1i*x0(i,2), x+1i*y]);
EtaV = 0; 
for i=1:N
    EtaV = EtaV + U(i)*subsij(C,i)*subsij0(normalization,i);
    EtaV = EtaV + U(i+N)*( subsij(Cx,i)*subsij0(normalization,i) + subsij(C,i)*subsij0(dx_norm,i));
    EtaV = EtaV + U(i+2*N)*( subsij(Cy,i)*subsij0(normalization,i)+ subsij(C,i)*subsij0(dy_norm,i) );
end
EtaV = symfun(EtaV,[x y]);

end