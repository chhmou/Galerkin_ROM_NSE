remap=0 * (1:NVU/2);



elnode2=elnode;
elnode2(:,4:6) = elnode2(:,4:6) + max(max(elnode(:,1:3)));

count=1;
for i=1:size(elnode2,1)
    for j=1:6
        if remap(elnode2(i,j))==0
            % new node
            remap(elnode2(i,j))=count;
            count=count+1;
        end
    end
end

remap2 = zeros(NVU,1);

remap2(1:2:end) = remap';
remap2(2:2:end) = remap'+NVU/2;

A = Acoeff(remap2,remap2);
A2= Acoeff2(remap2,remap2);

% reorder=[1:2:NVU,2:2:NVU];
% A = A(reorder,reorder);
% A2=A2(reorder,reorder);

spy(A)
figure
spy(A2)
