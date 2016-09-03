%% (0,1,1,0)-controlled Kingman f-seqs for n=5
f1=[1,0,0,1; 2,0,1,0; 3,1,0,0; 5,0,0,0;];
f2=[0,1,1,0; 2,0,1,0; 3,1,0,0; 5,0,0,0;];
f3=[0,1,1,0; 1,2,0,0; 3,1,0,0; 5,0,0,0;];

f=f3
fColswithPositiveSum=sum(f)>0;
fstar=f(:,fColswithPositiveSum);
% [U,W,V] = svd(fstar)
% w = 1 ./ diag(W,0)
% ww = zeros(size(w))
% Stablews=w>0.0001
% ww(Stablews)=w(Stablews)
% WW = diag(ww)
% Finv = V * WW * U'
f=fstar
t = [1,2,3,4]
Lambdas = 1 ./ [2*1/2, 3*2/2, 4*3/2, 5*4/2]
t = -(1 ./ Lambdas) .* log(rand(1,4))
l = t * f
f' * t'
Finv = pinv(f')
tEst = Finv * l'
(t-tEst') ./ t
 t*f
 tEst'*f


