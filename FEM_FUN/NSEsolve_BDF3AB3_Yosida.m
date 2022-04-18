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
        

        RHSvec(localUnk) = RHSvec(localUnk) + rhsf + rhscvc;
                
        
end

Acoeff =  [11/6/dt * MassMatrix + nu*StiffnessMatrix + gamma * GradDivMatrix,  PressureMatrix; PressureMatrix', spalloc(NPU,NPU,100*NPU) ];
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

% % Dirichlet pressure node
Acoeff(:,end)=1e-30 * Acoeff(:,end);
Acoeff(end,:)=1e-30 * Acoeff(end,:);
Acoeff(end,end)=1;
RHSvec(end)=0;

Acoeff = Acoeff .* (abs(Acoeff)>1e-16);

% Solve the linear system using Chorin/Teman splitting

% First, equip Mass matrix with boundary conditions
M = 11/6/dt * MassMatrix + 0*gamma*GradDivMatrix;
M(:,bdrydof) = 1e-30*M(:,bdrydof);
M(bdrydof,:) = 1e-30*M(bdrydof,:);
for ii=1:size(bdrydof,1)
        M(bdrydof(ii),bdrydof(ii))=1;
end
M = M .* (abs(M)>1e-16);

% We are using an inexact LU factorization
% U = [I M^-1 ; 0 I] and L = [A 0;Bt -Bt M^-1 B]

f = RHSvec(1:NVU);
g = RHSvec(NVU+1:end-1);
A = Acoeff(1:NVU,1:NVU);
B = Acoeff(1:NVU,NVU+1:end-1);
S = -B'*inv(M)*B;

tic 


% Solve Ly=b
y1 = A\f;
y2 = S \ (g - B'*y1);
% Solve Ux = b
GlobalP = [y2;0];
GlobalV = y1 - ( A \ (B*y2));

% Yosida = [ A, B ; B', B' * (inv(A) - inv(M))*B];
% Soln = Yosida \ RHSvec(1:end-1);
% 
% GlobalV2 = Soln(1:NVU);
% GlobalP2 = [Soln(NVU+1:end);0];
% 
% norm(GlobalV2 - GlobalV)

display(['Solve took ' num2str(toc) ' seconds'])

 
   
