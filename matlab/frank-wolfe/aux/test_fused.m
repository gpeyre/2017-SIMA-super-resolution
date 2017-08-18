N = 128;
D = spdiags([ones(N,1) -ones(N,1)],[0 -1],N,N); D = D(2:end,:);

rho = 20;

t = linspace(0,1,N)';
g = t - cos(6*pi*t)+.1;

cvx_solver SeDuMi % sdpt3 %  %
cvx_begin sdp quiet
cvx_precision high;
variable x(N,1);
norm(D*x,1) + rho*norm(x,1)/N <= 1;
maximize( sum(g.*x) )
cvx_end

clf; hold on;
plot(t, g, 'r:');
plot(t, .9*x/max(abs(x)), 'b');
axis tight;
box on;

%%
% Test for the proximal operator of |x|_1 s.t. sum(x)=0

proxL1 = @(x,lambda)sign(x).* max(0, abs(x) - lambda );

N = 10;
x = randn(N,1);
lambda = .5;
plot(proxL1(x,lambda), '.');

U = [];
for a=linspace(-10,10,1000)
    U(end+1) = sum( proxL1(x+a,lambda) );    
end


