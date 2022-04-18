%% Extrapolated BDF2 NSE solver
tic

i1 = zeros (100*Sdim, 1); 
j1 = zeros(100*Sdim, 1); 
mat1 = zeros(100*Sdim,1);
number1=0;

RHSvec = zeros(NVU+NPU,1) ;
RHSvec(1:NVU) = RHSvec(1:NVU) + 1/dt * MassMatrix * (2*PastV(:,end) - 0.5*PastV(:,end-1));

for itrg=1:size(elnode,1) 

        Vstart = [2*(elnode(itrg,1:3) - 1) + 1 , 2*(elnode(itrg,4:6) + nVert - 1) + 1];
        GlTrgVe = reshape([Vstart ; Vstart+1],[12,1]) ;            
        localUnk=GlTrgVe;
                
        % ( u^k-1 dot grad u^k,v)
        matcvc1 = inner_prod_ten1_VDiv_ten1(itrg,'quad_75','velFunBDF2', vel_bas_type, vel_bas_type) ;
         
        % RHS term 4: (f(t^n+1,v) 
        rhsf = inner_prod_ten1_Vec(itrg,'quad_75','rhsfun_bdf', vel_bas_type) ;
        
        RHSvec(localUnk) = RHSvec(localUnk) + rhsf ;
                
        for count=1:size(GlTrgVe,1)
           for count1=1:size(GlTrgVe,1)
                number1=number1+1;
                if(number1>length(i1)) 
                    error('init vecs bigger')
                end
                i1(number1) = GlTrgVe(count);
                j1(number1) = GlTrgVe(count1);
                mat1(number1) = matcvc1(count, count1);
           end
        end
end

NLMatrix = sparse(i1(1:number1), j1(1:number1), mat1(1:number1), NVU, NVU) ;
Acoeff =  [1.5/dt * MassMatrix + nu*StiffnessMatrix + gamma * GradDivMatrix + NLMatrix,  PressureMatrix; PressureMatrix', spalloc(NPU,NPU,100*NPU) ];
    
display(['Assembly took ' num2str(toc) ' seconds' ]) 


% Apply Dirichlet boundary conditions for velocity
bdrydirvals = getBdryVals(bdrydof,t);
RHSvec = RHSvec - Acoeff(:,bdrydof) * bdrydirvals;
Acoeff(:,bdrydof) = 1e-30*Acoeff(:,bdrydof);
Acoeff(bdrydof,:) = 1e-30*Acoeff(bdrydof,:);
for ii=1:size(bdrydof,1)
        Acoeff(bdrydof(ii),bdrydof(ii))=1;
end
RHSvec(bdrydof)=bdrydirvals;

% % Dirichlet pressure node
Acoeff(:,end)=1e-30 * Acoeff(:,end);
Acoeff(end,:)=1e-30 * Acoeff(end,:);
Acoeff(end,end)=1;
RHSvec(end)=0;

Acoeff = Acoeff .* (abs(Acoeff)>1e-16);



% Solve using Yosida splitting with pressure correction

% First, equip Mass matrix with boundary conditions
M = 1.5/dt * MassMatrix + gamma*GradDivMatrix + nu*GradDivMatrix;
M(:,bdrydof) = 1e-30*M(:,bdrydof);
M(bdrydof,:) = 1e-30*M(bdrydof,:);
for ii=1:size(bdrydof,1)
         M(bdrydof(ii),bdrydof(ii))=1;
end
M = M .* (abs(M)>1e-16);
N = Acoeff(1:NVU,1:NVU) - M;

f = RHSvec(1:NVU);
g = RHSvec(NVU+1:end-1);
A = Acoeff(1:NVU,1:NVU);
B = Acoeff(1:NVU,NVU+1:end-1);


 


% % i=1
% y1 = A\f;
% tic
% if gamma>0
%     [y2,FLAG,RELRES,ITER]= bicgstab(@(x)schurYosida2(x,B,M,N),-(g-B'*y1),1e-10,500,gamma*PressureMassMatrix(1:end-1,1:end-1) );
% else
%     [y2,FLAG,RELRES,ITER]= bicgstab(@(x)schurYosida2(x,B,M,N),-(g-B'*y1),1e-10,1000 );
% end
% display(['Flag = ' num2str(FLAG) ', iter = ' num2str(ITER) ', time=' num2str(toc) ])
% if FLAG>0
%       error('Schur complement solve failed')
% end
% GlobalP = [y2;0];
% GlobalV = y1 - ( A \ (B*y2));

% i=0
y1 = A\f;
tic
if gamma>0
    [y2,FLAG,RELRES,ITER]= pcg(@(x)schurYosida(x,B,M),-(g-B'*y1),1e-13,500,gamma*PressureMassMatrix(1:end-1,1:end-1) );
else
    [y2,FLAG,RELRES,ITER]= pcg(@(x)schurYosida(x,B,M),-(g-B'*y1),1e-10,1000 );
end
display(['Flag = ' num2str(FLAG) ', iter = ' num2str(ITER) ', time=' num2str(toc) ])
if FLAG>0
      error('Schur complement solve failed')
end
 GlobalP = [y2;0];
 GlobalV = y1 - ( A \ (B*y2));


% 
% % i=1 ???  (from quarteroni book)
% rhs1 = M \ (B*y2);
% rhs2 = B' * (M \ (A * rhs1));
% tic
% if gamma>0
%     [p,FLAG,RELRES,ITER]= pcg(@(x)schurYosida(x,B,M),rhs2,1e-13,500,PressureMassMatrix(1:end-1,1:end-1) );
% else
%     [p,FLAG,RELRES,ITER]= pcg(@(x)schurYosida(x,B,M),rhs2,1e-13,500 );
% end
% display(['Flag = ' num2str(FLAG) ', iter = ' num2str(ITER) ', time=' num2str(toc) ])
% if FLAG>0
%     error('Schur complement solve failed')
% end
% GlobalP = [p;0];
% GlobalV = y1 - ( A \ (B*p));
% % 


