
%UVSnap = Snapshots;

cutoff = eps;

% SnapD = Y^TM^hY/M
SnapD =UVSnap' * MassMatrix * UVSnap/(size(UVSnap,2));

[V,D] = eig(SnapD);

% find the cutoff
d = abs(diag(D));
e = d > d(1)*cutoff;
N = sum(e)

% do the svd
AlphaR = V(:,1:N);
D = diag(D);
LambdaR = D(1:N);

if averaging==1
PhiR0 = Uave;
end
  PhiR = zeros(NVU,N);
    for i=1:N
        x = UVSnap*AlphaR(:,i);
        PhiR(:,i) = x;
        denom = sqrt( x' * MassMatrix * x);
        PhiR(:,i) = PhiR(:,i) / denom;
    end


