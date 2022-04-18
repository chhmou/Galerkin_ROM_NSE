%% Extrapolated BDF2 NSE solver
tic

i1 = zeros (100*Sdim, 1); 
j1 = zeros(100*Sdim, 1); 
mat1 = zeros(100*Sdim,1);
number1=0;

RHSvec = zeros(NVU+NPU,1) ;
RHSvec(1:NVU) = RHSvec(1:NVU) + 1/dt * 0.5 * MassMatrix * (AllV(:,end) - AllV(:,end-1)) ...
                                - PressureMatrix*GlobalP ...
                                - nu*StiffnessMatrix*AllV(:,end);

for itrg=1:size(elnode,1) 

        Vstart = [2*(elnode(itrg,1:3) - 1) + 1 , 2*(elnode(itrg,4:6) + nVert - 1) + 1];
        GlTrgVe = reshape([Vstart ; Vstart+1],[12,1]) ;            
        localUnk=GlTrgVe;
                
        % ( u^k-1 dot grad u^k,v)
        matcvc1 = inner_prod_ten1_VDiv_ten1(itrg,'quad_75','velFunBDF2', vel_bas_type, vel_bas_type) ;
         
        % RHS term 4: (f(t^n+1,v) 
        rhsf = inner_prod_ten1_Vec(itrg,'quad_75','rhsfun_bdf', vel_bas_type) ;
        rhscvc = -inner_prod_ten1_Vec(itrg,'quad_75','velGradvelBDF2', vel_bas_type) ;
        
        
        RHSvec(localUnk) = RHSvec(localUnk) + rhsf + rhscvc;
                
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
bdrydirvals = getBdryValsUpdates(bdrydof,t);
RHSvec = RHSvec - Acoeff(:,bdrydof) * bdrydirvals;
Acoeff(:,bdrydof) = 1e-30*Acoeff(:,bdrydof);
Acoeff(bdrydof,:) = 1e-30*Acoeff(bdrydof,:);
for ii=1:size(bdrydof,1)
        Acoeff(bdrydof(ii),bdrydof(ii))=1;
end
RHSvec(bdrydof)=bdrydirvals;

% % Dirichlet pressure node
Acoeff(:,NVU+1)=1e-30 * Acoeff(:,NVU+1);
Acoeff(NVU+1,:)=1e-30 * Acoeff(NVU+1,:);
Acoeff(NVU+1,NVU+1)=1;
RHSvec(NVU+1)=0;

Acoeff = Acoeff .* (abs(Acoeff)>1e-16);






% % First, equip Mass matrix with boundary conditions
sig=1.5/dt;
M = gamma*GradDivMatrix + nu*StiffnessMatrix + sig*MassMatrix;
%M = sig*MassMatrix;
M(:,bdrydof) = 1e-30*M(:,bdrydof);
M(bdrydof,:) = 1e-30*M(bdrydof,:);
for ii=1:size(bdrydof,1)
         M(bdrydof(ii),bdrydof(ii))=1;
end
M = M .* (abs(M)>1e-16);

f = RHSvec(1:NVU);
g = RHSvec(NVU+1:end);
A = Acoeff(1:NVU,1:NVU);
B = Acoeff(1:NVU,NVU+1:end);

 

y1 = A\f;
tic
%if gamma==0
%    [y2,FLAG,RELRES,ITER]= pcg(@(x)schurYosida(x,B,M),-(g-B'*y1),1e-10,1000 );
%else
    [y2,FLAG,RELRES,ITER]= pcg(@(x)schurYosida(x,B,M),-(g-B'*y1),1e-10,1000,PressureMassMatrix );
%end
display(['Flag = ' num2str(FLAG) ', iter = ' num2str(ITER) ', time=' num2str(toc) ])
if FLAG>0
      error('Schur complement solve failed')
end
% GlobalP = [y2;0];
 GlobalPupdate = y2;
 GlobalVupdate = y1 - ( A \ (B*y2));

%  tic
%  soln = Acoeff \ RHSvec;
%  toc
% GlobalVupdate = soln(1:NVU);
% GlobalPupdate = soln(NVU+1:end);

% adjust the update
GlobalPupdate = GlobalPupdate - GlobalPupdate(1);

GlobalV = GlobalV + GlobalVupdate;
GlobalP = GlobalP + GlobalPupdate;
