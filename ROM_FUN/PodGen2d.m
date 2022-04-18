
function [PodU2dAll] = PodGen2d(UVSnap,Uave,MassMatrix,ModeD,averaging)
NVU = size(UVSnap,1);

cutoff = eps;
% if averaging==1
%     Uave = sum(UVSnap')/size(UVSnap,2);
%     Uave = Uave';
%     for i=1:size(UVSnap,2)
%          UVSnap(:,i) = UVSnap(:,i) - Uave;
%     end
% end

%U =U/sqrt(size(U,2));

%UVSnap = UVSnap/sqrt(size(UVSnap,2));
M =UVSnap' * MassMatrix * UVSnap/(size(UVSnap,2));

[V,D] = eig(M);

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


if averaging==1
    PodU2dAll = zeros(length(Uave),ModeD);
    PodU2dAll(:,2:ModeD) = PhiR(:,1:ModeD-1);
    PodU2dAll(:,1) = Uave;
else
    PodU2dAll = PhiR(:,1:ModeD);
end
