%% Extrapolated BDF2 NSE solver
tic

i1 = zeros (100*Sdim, 1); 
j1 = zeros(100*Sdim, 1); 
mat1 = zeros(100*Sdim,1);
number1=0;

RHSvec = zeros(NVU+NPU,1) ;
RHSvec(1:NVU) = RHSvec(1:NVU) + 1/dt * MassMatrix * (2*AllV(:,end) - 0.5*AllV(:,end-1)); %% BDF2 discretization move to the right

for itrg=1:size(elnode,1) 

        Vstart = [2*(elnode(itrg,1:3) - 1) + 1 , 2*(elnode(itrg,4:6) + nVert - 1) + 1];
        GlTrgVe = reshape([Vstart ; Vstart+1],[12,1]) ;            
        localUnk=GlTrgVe;
                
        ten1a=basisValues(:,:,triag_no);
        Gradtrue = gradBasisValues(:,:,:,triag_no);
        
        % ( u^k-1 dot grad u^k,v)
        matcvc1 = inner_prod_ten1_VDiv_ten1(itrg,'quad_75','velFunBDF2', vel_bas_type, vel_bas_type) ;
         
%        % RHS term 4: (f(t^n+1,v) 
%        rhsf = inner_prod_ten1_Vec(itrg,'quad_75','rhsfun_bdf', vel_bas_type) ;
        
%        RHSvec(localUnk) = RHSvec(localUnk) + rhsf ;
                
        if abs(t-dt)<1e-8
        
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
        else
            NLMatrix(localUnk,localUnk)=NLMatrix(localUnk,localUnk)+matcvc1;
        end
end
if abs(t-dt)<1e-8
    NLMatrix = sparse(i1(1:number1), j1(1:number1), mat1(1:number1), NVU, NVU) ;
end
Acoeff =  [1.5/dt * MassMatrix + nu*StiffnessMatrix + gamma * GradDivMatrix + NLMatrix,  PressureMatrix; PressureMatrix', spalloc(NPU,NPU,100*NPU) ];

NLMatrix = 1e-30*NLMatrix;

display(['Assembly took ' num2str(toc) ' seconds' ]) 


% Apply Dirichlet boundary conditions for velocity
bdrydirvals = GlobalV(bdrydof);
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

Acoeff = Acoeff .* (abs(Acoeff)>1e-16); %%MOU%% why write this (abs(Acoeff)>1e-16)

tic
soln = Acoeff \ RHSvec;
toc 
GlobalV = soln(1:NVU);
GlobalP = soln(NVU+1:end);
