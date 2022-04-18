%% Extrapolated BDF3 NSE solver (Karniadakis et al 1991)

tic

RHSvec = zeros(NVU+NPU,1) ;
RHSvec(1:NVU) = RHSvec(1:NVU) + 1/dt * MassMatrix * (3*PastV(:,end) - 1.5*PastV(:,end-1) + 1/3*PastV(:,end-2) );
%RHSvec(1:NVU) = RHSvec(1:NVU) + 1/dt * MassMatrix * (2*PastV(:,end) - 0.5*PastV(:,end-1));

for itrg=1:size(elnode,1) 

        Vstart = [2*(elnode(itrg,1:3) - 1) + 1 , 2*(elnode(itrg,4:6) + nVert - 1) + 1];
        GlTrgVe = reshape([Vstart ; Vstart+1],[12,1]) ;            
        localUnk=GlTrgVe;
                
        % RHS term 1: - 3 * ( u^n dot grad u_n,v) 
        rhscvc = -3 * inner_prod_ten1_Vec(itrg,'quad_75','velPGradvelP', vel_bas_type) ;

        % RHS term 2: + 3 * ( u^n-1 dot grad u_n-1,v) 
        rhscvc = rhscvc + 3*inner_prod_ten1_Vec(itrg,'quad_75','velPPGradvelPP', vel_bas_type) ;

        % RHS term 3: -  ( u^n-2 dot grad u_n-2,v) 
        rhscvc = rhscvc - inner_prod_ten1_Vec(itrg,'quad_75','velPPPGradvelPPP', vel_bas_type) ;
        
        % RHS term 4: (f(t^n+1,v) 
        rhsf = inner_prod_ten1_Vec(itrg,'quad_75','rhsfun_bdf', vel_bas_type) ;
        

        RHSvec(localUnk) = RHSvec(localUnk) + rhsf;
                
        
end

Acoeff =  [11/6/dt * MassMatrix + nu*StiffnessMatrix ,  PressureMatrix; PressureMatrix', spalloc(NPU,NPU,100*NPU) ];
%Acoeff =  [1.5/dt * MassMatrix + nu*StiffnessMatrix + gamma * GradDivMatrix,  PressureMatrix; PressureMatrix', spalloc(NPU,NPU,100*NPU) ];
    
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

% Dirichlet pressure node
Acoeff(:,end)=1e-30 * Acoeff(:,end);
Acoeff(end,:)=1e-30 * Acoeff(end,:);
Acoeff(end,end)=1;
RHSvec(end)=0;

Acoeff = Acoeff .* (abs(Acoeff)>1e-16);

% Solve the linear system   
tic
Ukn = Acoeff \ RHSvec ;
display(['Solve took ' num2str(toc) ' seconds'])

 GlobalV = Ukn(1:NVU) ;
 GlobalP = Ukn( NVU+1:end ) ;
 
 % adjust GlobalP -> 1st entry should be sin(1+t);
 %GlobalP = GlobalP - GlobalP(1) + sin(1+t);
 
   
