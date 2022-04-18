%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%

%#function CtsQuad CtsLin DCtsLin quad_75 velFunPrevTime velFun
%#function GradvelFunPrevTimedef GradvelFun GradvelFunPrevTime
%#function velGradvel veltimesvelFun velPGradvelP velPtimesvelPFun
%#function onefun Identity OneScalFun extrapvelGradvelPrev extrapvelFun

%%%%%%% In this routine we investigate the use of the Scott-Vogelius elements
%%% for the Navier Stokes problem  --- 09/07/09
%
%
%
%clear

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  GlobalS  GlobalG nu PastV problem
global dimTvel  dimTpre  dimTstr  dimTGrv 
global vel_bas_type  pre_bas_type  str_bas_type  Grv_bas_type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Define the problem and the solution method
%%

% We begin by specifying the type of approximation elements
% we are using for the velocity, pressure and the overall solution method.

%
vel_bas_type = 'CtsQuad' ;
pre_bas_type = 'DCtsLin' ;

gd_parameter = 0;

%problem = 'Cylinder'
problem = 'Step'


if strcmp(problem,'Cylinder')
    nu=0.001;
    t=0;
    T=8;
    dt=0.01;
elseif strcmp(problem,'Step')
    nu=1/600;
    t=0;
    T=40;
    dt=0.01;
else
    error('Please specify a valid problem to solve')
end


numTimeSteps = round(T/dt);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Set up the Geometry and the Boundary/Initial conditions.

% We assume that the subroutine will pass back 
% nodeco, elnode, and bdynde.
% nodeco will contain the x and y coordinate of the nodes
% elnode will contain the definition of the triangulation
% bdynde will contain the definition of the boundary in a
% counter clockwise direction. 

% Enter the grid file
meshlevel=35;

% with barycenter refinement
GenSVgrd('Stepmesh3b.msh');
[nodeco, elnode, bdynde] = FFEM2ErvGrd2('SVgrid.msh') ;

%% w/out barycenter refinement
%[nodeco, elnode, bdynde] = FFEM2ErvGrd2('CylinderMesh1b.msh') ;


%% Set up the Dirichlet B.C.s for the velocity. 
%% Remark -- We also set up a "stress" bc which is not used -- needed for midEdgegen
nbdy = size(bdynde,1) ;
bdynde = [bdynde , zeros(nbdy , 1), -1*ones(nbdy , 2) ] ;
[nodeco, elnode, bdynde, bdyedge, nVert, nedge] =  midEDGEgen(nodeco, elnode, bdynde) ;

% If no partitioning of the region is desired, use the following
% commands.
 ntri = size(elnode,1) ;
 DomPoint(1) = 1 ;
 TrgMeNxt = 2:1:ntri ;
 TrgMeNxt(ntri) = -1 ;
%[DomPoint, TrgMeNxt] = domsubstr(nodeco, elnode) ;




%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Set up and initialize the Global Solution vectors.
%
% Note associated with each unknown in the global solution vector
% is a indicator which denotes if the unknowns represents a 
% Dirichlet boundary condition.


nDom = size(DomPoint,2) ;

if strcmp(vel_bas_type, 'CtsQuad') == 1 
   GlobalV = zeros((nVert + nedge)*2, 2) ;
   GlobalV(:,2) = -1 ;

   for ibdy = 1 : size(bdynde,1)
      if ( bdynde(ibdy,2) ~= -1 )
         GlobalV( 2*bdynde(ibdy,1) - 1 , 2 ) = bdynde(ibdy,2) ;
         GlobalV( 2*bdynde(ibdy,1) , 2 ) = bdynde(ibdy,2) ;
      end
   end
   
   for ibdye = 1 : size(bdyedge,1)
      if ( bdyedge(ibdye,2) ~= -1 )
         GlobalV( 2* nVert + 2*bdyedge(ibdye,1) - 1 , 2 ) = bdyedge(ibdye,2) ;
         GlobalV( 2* nVert + 2*bdyedge(ibdye,1) , 2 ) = bdyedge(ibdye,2) ;     
      end
   end
   
   if strcmp(problem,'Cylinder')
       % For the cylinder problem, something is wrong with setting bc via
       % bdyedge - it misses points.  Hence we will do this directly
       for ii=1:size(nodeco,1)
           xypt = nodeco(ii,:);
           if abs( xypt(1)*(xypt(1)-2.2) )<1e-10 | abs( xypt(2)*(xypt(2) - .41) )<1e-10 | ...
                   sqrt( (xypt(1) - .2)^2 + (xypt(2)-.2)^2 ) < 0.05001
               GlobalV(2*ii-1 : 2*ii ,2)=0;
           end
       end
   end
           
   
   % Set the boundary/ initial velocity. 
   tempV = initVel(nodeco) ;
   GlobalV(1:2:(nVert + nedge)*2,1) = tempV(1,:).' ;
   GlobalV(2:2:(nVert + nedge)*2,1) = tempV(2,:).' ;
   

end


if strcmp(pre_bas_type, 'CtsLin') == 1 
   GlobalP = zeros(nVert,2) ;
   GlobalP(:,2) = -1 ;
   for ibdy = 1 : size(bdynde,1)
      if ( bdynde(ibdy,3) == 0 )
         GlobalP( bdynde(ibdy,1) , 2 ) = 0 ;
      end
   end
   
elseif strcmp(pre_bas_type, 'DCtsLin') == 1 
   ntri = size(elnode,1) ;
   GlobalP = zeros(3*ntri,2) ;
   GlobalP(:,2) = -1 ;
   mapN2GP = reshape(elnode(:,1:3).',[3*ntri 1]) ;
   for ibdy = 1 : size(bdynde,1)
      if ( bdynde(ibdy,3) ~= -1 )
         for inde = 1 : 3*ntri
            if ( bdynde(ibdy,1) == mapN2GP(inde) )
               GlobalP( inde , 2 ) = bdynde(ibdy,3) ;
            end
         end
      end
   end
   
end

NVU = size(GlobalV,1)
NPU = size(GlobalP,1)
NVU+NPU
% Get boundary info for ease in later boundary condition applications
% Assumes all dof of a node are dirichlet if one is.
bdrydof=[];
for ii=1:size(GlobalV,1)
    if GlobalV(ii,2)==0
        bdrydof = [bdrydof; ii];
    end
end




%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Now to solve the problem !  

%% create the matrices that will be involved in each timestep
createConstantMatrices2;

%% Now go through the time steps and solve!

% % Put initial condition in Vh
% Acoeff = [MassMatrix, PressureMatrix; PressureMatrix', spalloc(NPU,NPU,10)];
% RHSvec = [MassMatrix*GlobalV(:,1) ; zeros(NPU,1)];
% 
% %% Apply Dirichlet boundary conditions for velocity
% bdrydirvals = getBdryVals(bdrydof,0);
% RHSvec = RHSvec - Acoeff(:,bdrydof) * bdrydirvals;
% Acoeff(:,bdrydof) = 1e-30*Acoeff(:,bdrydof);
% Acoeff(bdrydof,:) = 1e-30*Acoeff(bdrydof,:);
% for ii=1:size(bdrydof,1)
%     Acoeff(bdrydof(ii),bdrydof(ii))=1;
% end
% RHSvec(bdrydof)=bdrydirvals;
% 
% 
% %% Replace the last pressure equation with mean pressure = 0
% Acoeff(Sdim, :) = pconst ;
% RHSvec(Sdim) = 0 ;        
% 
% 
% xx = Acoeff \ RHSvec;
% GlobalV(:,1) = xx(1:NVU);

% Try SV u0 as u0 for TH
 if meshlevel==35 & strcmp(problem,'Step')
     load StepSV_u0_mesh3b
 end

%% PastV will hold the three previous time steps
PastV = [GlobalV(:,1),GlobalV(:,1),GlobalV(:,1)];
PastP = [GlobalP(:,1),GlobalP(:,1),GlobalP(:,1)];


CalcDivErr


for jj=1:numTimeSteps
    
    t = jj*dt
    
%    NSEsolve_nonlinear_SV
    NSEsolve_CNLE_SV
    
    
    % Uses GlobalV data to calculate ||div u_h||
    CalcDivErr
    divL2
    divL2err(jj)=divL2;

    if mod(jj,500)==0 || jj==10 || jj==100
    
        filename = ['Step_NSE_SV_level' num2str(meshlevel) '_g' num2str(gd_parameter) '_timestep' num2str(jj)]
        save(filename,'GlobalV','nodeco','GlobalP')
    end
    
    PastV = [PastV(:,2:3),GlobalV(:,1)];
end

divL2err

      
