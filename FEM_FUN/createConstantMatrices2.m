% Creates the matrices for (u,v), (D(u),D(v)), (p,div v) and (div u,v) and (grad u,grad v)

% This is a script, so I have access to all info, like boundary info, mesh info, etc

%% Guess an allocation for the matrices
% The size of the Stiffness and Mass are only to be the upper left corner -
% they will be used for filtering as well as NSE solve.







tic;

i1 = zeros (50*Sdim, 1); 
j1 = zeros(50*Sdim, 1); 
i2 = zeros (50*Sdim, 1); 
j2 = zeros(50*Sdim, 1); 
i3 = zeros (50*Sdim, 1); 
j3 = zeros(50*Sdim, 1); 


pm1 = zeros(50*Sdim,1);
p1 = zeros(50*Sdim, 1);
c1= zeros(50*Sdim, 1);
s1 = zeros(50*Sdim, 1);
s2 = zeros(50*Sdim, 1); 
m1 = zeros(50*Sdim, 1);
g1 = zeros(50*Sdim,1);

number1=0;
number2=0;
number3=0;

% nUkn = 0 ;
% mapGV2Uk = zeros(size(GlobalV,1),1).' ;
% mapUk2GV = [ ] ;
% mapGP2Uk = zeros(size(GlobalP,1),1).' ;
% mapUk2GP = [ ] ;

pconst = zeros(1, NVU + NPU) ;

for itrg = 1:size(elnode,1)
            if mod(itrg,100)==1
%                display(num2str(itrg))
            end
      
          Alocal = [ ] ;
          RHSloc = [ ] ;
          localUnk = [ ] ;
          ploc = [ ] ;
          % Set up unknown solution vector mapping to global velocity vector
          if strcmp(vel_bas_type, 'CtsQuad') == 1
                Vstart = [2*(elnode(itrg,1:3) - 1) + 1 , 2*(elnode(itrg,4:6) + nVert - 1) + 1];
                GlTrgVe = reshape([Vstart ; Vstart+1],[12,1]) ;

          end


          % Set up unknown solution vector mapping to global pressure vector
          if strcmp(pre_bas_type, 'CtsLin') == 1
                Pstart = elnode(itrg, 1:3) ;
                GlTrgPr = Pstart.' ;

          elseif strcmp(pre_bas_type, 'DCtsLin') == 1
             GlTrgPr = [3*itrg - 2 ; 3*itrg - 1 ; 3*itrg ] ;

          end

          localUnk = [ GlTrgVe ; NVU + GlTrgPr ] ;

          %%%%%%%%%%%%%%%%
          % Evaluate the integrals.

          %% Assemble the pieces for the MOMENTUM EQUATION

          %% Local Stiffness (in deformation form)
          matA11_a = inner_prod_Grad_ten1_Grad_ten1(itrg, 'quad_75', ...
             'onefun', vel_bas_type, vel_bas_type) ;
          matA11_b = inner_prod_Grad_ten1T_Grad_ten1(itrg, 'quad_75', ...
             'onefun', vel_bas_type, vel_bas_type) ;
          LStiff = 0.5 * (matA11_a + matA11_b);

          % Local Mass
          LMass = inner_prod_ten1_ten1_SqM(itrg,'quad_75', ...
              'Identity', vel_bas_type, vel_bas_type) ;

          % Local Pressure Mass
          LPresMass = inner_prod_ten0_ten0(itrg,'quad_75',pre_bas_type, pre_bas_type) ;

          % Local Grad div
          LGD = inner_prod_divten1_divten1(itrg,'quad_75', ...
              'OneScalFun', vel_bas_type, vel_bas_type) ;

          % Pressure (and cons mass if you transpose)
          LPres = -1.0* inner_prod_Div_ten1_ten0(itrg,'quad_75', ...
             'OneScalFun', vel_bas_type, pre_bas_type) ;     

%           %%% Sept 15: Calculation to impose mean pressure zero constraint
%           ploc = inner_prod_ten0(itrg, 'quad_75', ...
%              'OneScalFun', pre_bas_type) ;
% 
% 
% 
%           % Store these calcuations.
%           pconst(NVU + [GlTrgPr]) = pconst(NVU + [GlTrgPr]) + ploc.' ; 


          for count=1:size(GlTrgVe,1)
           for count1=1:size(GlTrgPr,1)
                number1=number1+1;
                if(number1>length(i1)) 
                    error('init i1 to be bigger')
                    i1 (2*size(i1,1))=0; 
                    j1 (2*size(i1,1))=0; 
                    p1 (2*size(i1,1))=0;
                    c1 (2*size(i1,1))=0;
                end
                i1(number1) = GlTrgVe(count);
                j1(number1) = GlTrgPr(count1);
                p1(number1) = LPres(count, count1);
                c1(number1) = LPres(count, count1)';    
           end
          end

        for count=1:size(GlTrgVe,1)
           for count1=1:size(GlTrgVe,1)
                number2=number2+1;
                if(number2>length(i2)) 
                    error('init i2 to be bigger')
                    i2 (2*size(i2,1))=0;
                    j2 (2*size(j2,1))=0;
                    s1 (2*size(i2,1))=0;
                    s2 (2*size(i2,1))=0;
                    m1 (2*size(i2,1))=0;
                end
                i2(number2) = GlTrgVe(count);
                j2(number2) = GlTrgVe(count1);
                s1(number2) = LStiff(count, count1);
                m1(number2) = LMass(count, count1);
                s2(number2) = matA11_a(count, count1);
                g1(number2) = LGD(count,count1);

           end
        end

        for count=1:size(GlTrgPr,1)
          for count1=1:size(GlTrgPr,1)
                number3=number3+1;
                if(number3>length(i3)) 
                    error('init i3 to be bigger')
                end
                i3(number3) = GlTrgPr(count);
                j3(number3) = GlTrgPr(count1);
                pm1(number3) = LPresMass(count, count1);

           end
        end


   
end

PressureMatrix = sparse(i1(1:number1), j1(1:number1), p1(1:number1), NVU,NPU) ;
ConsMassMatrix = sparse(j1(1:number1), i1(1:number1), c1(1:number1), NPU,NVU);
StiffnessMatrix = sparse(i2(1:number2), j2(1:number2), s1(1:number2), NVU, NVU) ;
MassMatrix = sparse(i2(1:number2), j2(1:number2), m1(1:number2), NVU,NVU) ;
StiffnessMat2 = sparse(i2(1:number2), j2(1:number2), s2(1:number2), NVU, NVU) ;
GradDivMatrix = sparse(i2(1:number2), j2(1:number2), g1(1:number2),NVU,NVU);
PressureMassMatrix = sparse(i3(1:number3), j3(1:number3), pm1(1:number3),NPU,NPU);

%display(['Initial matrix assemblies took ' num2str(toc) ' seconds.'])

StiffnessMatrix = StiffnessMat2;

clear i1 i2 i3 j1 j2 j3 p1 c1 s1 m1 s2 g1 pm1

