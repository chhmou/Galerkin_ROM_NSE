%% Extrapolated BDF2 NSE solver
tic

i1 = zeros (100*Sdim, 1); 
j1 = zeros(100*Sdim, 1); 
mat1 = zeros(100*Sdim,1);
number1=0;

RHSvec = zeros(NVU,1) ;

for itrg=1:size(elnode,1) 

        Vstart = [2*(elnode(itrg,1:3) - 1) + 1 , 2*(elnode(itrg,4:6) + nVert - 1) + 1];
        GlTrgVe = reshape([Vstart ; Vstart+1],[12,1]) ;            
        localUnk=GlTrgVe;
                
        % ( U dot grad u^k,v)
        matcvc1 = inner_prod_ten1_VDiv_ten1(itrg,'quad_75','velFun', vel_bas_type, vel_bas_type) ;
         
        % RHS term 4: (f(t^n+1,v) 
        rhsf = inner_prod_ten1_Vec(itrg,'quad_75','rhsfun', vel_bas_type) ;
        
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
%Acoeff =  [nu*StiffnessMatrix + gamma * GradDivMatrix + NLMatrix,  PressureMatrix; PressureMatrix', sparse(NPU,NPU) ];
Acoeff =  [sig*MassMatrix + nu*StiffnessMatrix + gamma * GradDivMatrix + NLMatrix];
    
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

Acoeff = Acoeff .* (abs(Acoeff)>1e-16);

% %direct solve
% x = Acoeff \ RHSvec;
% GlobalV = x(1:NVU);

% iterated penalty method
% step 0

format shorte
table=[];
sumV = 0*GlobalV;
for i=1:2
    RHSvec2 = RHSvec - gamma*GradDivMatrix*sumV;
    RHSvec2(bdrydof)=bdrydirvals;
    GlobalV = Acoeff \ RHSvec2;
    sumV = sumV + GlobalV;
    % calculate error
    CalcErrDiv
    % calculate difference
    temp=GlobalV;
    temp5 = GlobalV;
    GlobalV = temp - GlobalVSV;
    CalcErrDiv3
    GlobalV = temp5;
    table = [table;H1norm,divL2,preErr3,GradVelL2Error]
    display('****************************')
end
    
    
    


