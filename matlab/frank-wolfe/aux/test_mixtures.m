rep = 'results/mixtures/';
if not(exist(rep))
    mkdir(rep);
end

% with wrong normalization: non-degenerated ==> BUG
H = @(x,z)1./sqrt(z) .* exp( -x.^2 .* z );
Hz = @(x,z)-(x.^2 + 1./(2*z)).*H(x,z);
Hxx = @(x,z)2*z.*(2*x.^2.*z - 1).*H(x,z);
Hzz = @(x,z)( 3/(4*z.^2) + x.^2./(2*z) + x.^4 ).*H(x,z);

% with correct normalization: degenerated
H = @(x,z)sqrt(z) .* exp( -x.^2 .* z );
Hz = @(x,z)(1./(2*z) - x.^2 ).*H(x,z);
Hxx = @(x,z)2*z.*(2*x.^2.*z - 1).*H(x,z);
Hzz = @(x,z)( x.^4 - x.^2./z - 1./(4*z.^2) ).*H(x,z);

% with no normalization: singular Gamma
H = @(x,z)exp( -x.^2 .* z );
Hz = @(x,z)-x.^2.*H(x,z);
Hxx = @(x,z)( 4*x.^2.*z.^2-2*z ).*H(x,z);
Hzz = @(x,z)x.^4.*H(x,z);


z = linspace(.01,2,512);
x = linspace(-3,3,513);
[Z,X] = meshgrid(z,x);


x0=0; z0=1/2;
Gamma = [H(0,2*z0), 0, Hz(0,2*z0);  0, Hxx(0,2*z0), 0; Hz(0,2*z0), 0, Hzz(0,2*z0) ];
U = Gamma \ [1;0;0]; % H(2) must be 0
eta = @(x,z)U(1)*H(x,z+1/2) + U(3)*Hz(x,z+1/2);

A = eta(X,Z)';
figure(1); 
clf; hold on;
imagesc(x,z,A);
colormap parula(256);
contour(x,z,A, linspace(-1,.999,20), 'k');
box on;
saveas(gcf, [rep 'mixtures-2d.png'], 'png');

figure(2); 
clf;
subplot(2,1,1);
plot(x, eta(x,1/2));
axis tight;
subplot(2,1,2);
plot(z, eta(0,z));
axis tight;
saveas(gcf, [rep 'mixtures-2d.eps'], 'epsc');


a = [1 0 -1/2; 0 -2 0; -1/2 0 3/4];