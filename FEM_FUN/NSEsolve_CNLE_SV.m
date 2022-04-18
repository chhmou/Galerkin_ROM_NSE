%% Extrapolated BDF3 NSE solver (Karniadakis et al 1991)

tic
%% Assemble and solve
i2 = zeros (200*Sdim, 1); 
j2 = zeros(200*Sdim, 1); 
uvmat = zeros(200*Sdim, 1);
number2=0;

RHSvec = zeros(NVU+NPU,1) ;
RHSvec(1:NVU) = RHSvec(1:NVU) + 1/dt * MassMatrix * (3*PastV(:,end) - 1.5*PastV(:,end-1) + 1/3*PastV(end-2) );

for itrg=1:size(elnode,1) 

        Vstart = [2*(elnode(itrg,1:3) - 1) + 1 , 2*(elnode(itrg,4:6) + nVert - 1) + 1];
        GlTrgVe = reshape([Vstart ; Vstart+1],[12,1]) ;            
        localUnk=GlTrgVe;
                
        % RHS term 1: - 3 * ( u^n dot grad u_n,v) 
        rhscvc = -3 * inner_prod_ten1_Vec(itrg,'quad_75','velPGradvelP', vel_bas_type) ;

        % RHS term 2: + 3 * ( u^n-1 dot grad u_n-1,v) 
        rhscvc = rhscvc + 3*inner_prod_ten1_Vec(itrg,'quad_75','velPPGradvelPP', vel_bas_type) ;

        % RHS term 1: -  ( u^n-2 dot grad u_n-2,v) 
        rhscvc = rhscvc - inner_prod_ten1_Vec(itrg,'quad_75','velPPPGradvelPPP', vel_bas_type) ;

        
        %% Assemble the local matrices
        RHSvec(localUnk) = RHSvec(localUnk) + rhscvc + rhsgd;
                
        
        for count=1:size(GlTrgVe,1)
           for count1=1:size(GlTrgVe,1)
                number2=number2+1;
                if(number2>length(i2)) 
                    error('Need to initialize vecs bigger in NSEsolve_nonlinear_SV')
                end
                i2(number2) = GlTrgVe(count);
                j2(number2) = GlTrgVe(count1);
                uvmat(number2) = matcvc(count, count1);
           end
        end

             
        itrg = TrgMeNxt(itrg) ;  
       end

    end
    UVMatrix = sparse(i2(1:number2), j2(1:number2), uvmat(1:number2), NVU,NVU) ;
    Acoeff =  [1/dt * MassMatrix + nu/2*StiffnessMatrix + gd_parameter/2 * GradDivMatrix + UVMatrix,  PressureMatrix; ConsMassMatrix, spalloc(NPU,NPU,100*NPU) ];
    
    display(['Assembly took ' num2str(toc) ' seconds' ]) 

    %% Apply Dirichlet boundary conditions for velocity
    bdrydirvals = getBdryVals(bdrydof,t);
    RHSvec = RHSvec - Acoeff(:,bdrydof) * bdrydirvals;
    Acoeff(:,bdrydof) = 1e-30*Acoeff(:,bdrydof);
    Acoeff(bdrydof,:) = 1e-30*Acoeff(bdrydof,:);
    for ii=1:size(bdrydof,1)
        Acoeff(bdrydof(ii),bdrydof(ii))=1;
    end
    RHSvec(bdrydof)=bdrydirvals;
    
    
    %% Replace the last pressure equation with mean pressure = 0
    Acoeff(Sdim, :) = pconst ;
    RHSvec(Sdim) = 0 ;         

    %% Solve the linear system   
    tic
    Ukn = Acoeff \ RHSvec ;
    display(['Solve took ' num2str(toc) ' seconds'])
    
    GlobalV(:, 1) = Ukn(1:NVU) ;
    GlobalP(:, 1) = Ukn( NVU+1:NVU+NPU ) ;
   
