%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%

%%%%%%% In this routine we investigate the use of the Scott-Vogelius elements
%%% for the Stokes problem  --- 09/07/09
%
%  This is the driver routine for a Stokes calculation
%
%
clear

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  GlobalS  GlobalG nu problem
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

nu=1/100;
hhh=64;


problem='cavity'


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
gridfile=['SQ' num2str(hhh) '.msh']
GenSVgrd(gridfile);
[nodeco, elnode, bdynde] = FFEM2ErvGrd2('SVgrid.msh') ;

%%% ignore bdynde
elnode2 = [nodeco(elnode(:,1),1), nodeco(elnode(:,1)), elnode];
elnode2=sortrows(elnode2);


%[nodeco, elnode, bdynde] = FFEM2ErvGrd2('SQ16.msh') ;


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
Sdim = NVU+NPU


% Get boundary info for ease in later boundary condition applications
% Assumes all dof of a node are dirichlet if one is.
bdrydof=[];
for ii=1:size(GlobalV,1)
    if GlobalV(ii,2)==0
        bdrydof = [bdrydof; ii];
    end
end


%% PastV will hold the three previous time steps
PastV = [GlobalV(:,1),GlobalV(:,1),GlobalV(:,1)];
PastP = [GlobalP(:,1),GlobalP(:,1),GlobalP(:,1)];

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Now to solve the problem !  

%% create the matrices that will be involved in each timestep

if hhh==16
    load constantmats_cavity_16    
elseif hhh==32
    load constantmats_cavity_32
elseif hhh==64
    load constantmats_cavity_64
elseif hhh==96
    load constantmats_cavity_96
else
    createConstantMatrices2;
    save constantmats_cavity_128 MassMatrix StiffnessMatrix PressureMatrix GradDivMatrix PressureMassMatrix
end

%load dc400soln128 GlobalV
NSEsolve_nonlinear_steady_SV_conv

CalcErrDiv



