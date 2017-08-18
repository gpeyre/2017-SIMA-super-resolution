a = 0.1;
op = load_operator('laplace-normalized',a);

x = linspace(0,3,1000);
x0 = 1;
lambda = .7; 


q = 6;
lambda = linspace(0,2,q);
lambda = [.3 .5  .6 .8 1 1.5 3];
q = length(lambda);
clf; hold on;
for i=1:q
    col = [(i-1)/(q-1), 0,1-(i-1)/(q-1)];
    plot(x, z(lambda(i)), 'color', col);
    lgd{i} = num2str(lambda(i)); 
end
axis tight; box on;
legend(lgd);

clf;
plot(x, z(.6));