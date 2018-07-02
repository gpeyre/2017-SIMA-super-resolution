%Run:
%>syms u1 u2 v1 v2;
%>[C,L,M] = GaussElim1([u1 u2; v1 v2],10);
%C is the de Boor basis
%L is the matrix L_{tz} from 2spikes_v2.pdf
%M is the matrix L_{tz}^{-1} V_{tz}

function [C,L,M] = GaussElim1(X,N)
% % N=4;
% % X= round(randn(N,2)*10);
% % N=2*N-1;

syms t u v a b

syms x y

Nodes =sym(X);% sym(round(X));
% Nodes = [0 0;  1 0; u 0];

Nodes = t*Nodes; % Nodes(1,:) = sym([0 0]);

% Nodes = [0 0;  1 0;u v; 0 1];
len_n = length(Nodes(:,1));

%Create the Taylor expansion matrix
M=[]; i = 1; colInfo=sym('A', [(N+2)*(N+1)/2  1]);

M0 = []; M1 = []; M2 = [];
for p = 1:len_n
    r0 = [];
    r11 = [];
    r12 = [];
    for n=0:N
        for j=0:n
            k= n-j;
            if p==1
                colInfo(i) = x^j*y^k; i=i+1;
            end
            r0 = [r0 taylorcoeff(Nodes(p,1), Nodes(p,2), j,k)];
            %             if p~=len_n
            r11 = [r11 taylorcoeff(Nodes(p,1), Nodes(p,2), j-1,k)];
            r12 = [r12 taylorcoeff(Nodes(p,1), Nodes(p,2), j,k-1)];
            %             end
        end
    end
    %     M = [M;r0;r11;r12];
    M0 = [M0;r0];
    M1 = [M1;r11];
    M2 = [M2;r12];
end
M = [M0;M1;M2];
M_taylor = M;



%Get the row reduced echelon form
M = rref(M_taylor);
M = simplify(M);

ll = size(M_taylor,2);
M2 = [M_taylor, sym(eye(3*len_n))];
M3 = rref(M2);
L=M3(:,ll+1:end);

C= M*colInfo;

%select the least degree homogeneous polynomials
for j=1:length(C)
    c_tmp = coeffs(C(j),t);
    C(j) = c_tmp(1);
end



end


%get taylor coefficient of order (n,m), centred at (0,0) in direction (u,v)
function coeff = taylorcoeff(u,v, n,m)
if n>=0
    un = u^n;
    fn = factorial(n);
else
    un = 0;  fn = 1;
end

if m>=0
    vm = v^m;  fm = factorial(m);
else
    vm = 0;  fm = 1;
end

coeff = un*vm/fn/fm;
end
