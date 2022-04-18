%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDF3 + algebraic splitting


clear
format long

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge GlobalVnew
global GlobalV  GlobalP  nu vel_bas_type  pre_bas_type problem sig
global GlobalPSV sumV t AllV




vel_bas_type = 'CtsQuad' ;  
pre_bas_type = 'DCtsLin' ;  % could be DCtsLin

nu=.1;
t=0;

T=1.0;
dt=T/4;
    
problem = 'YosidaTest2'

meshname = 'SquareMeshU16.msh'
GenSVgrd(meshname);
[nodeco, elnode, bdynde] = FFEM2ErvGrd2('SVgrid.msh') ;


nbdy = size(bdynde,1) ;
bdynde = [bdynde , zeros(nbdy , 1), -1*ones(nbdy , 2) ] ;
[nodeco, elnode, bdynde, bdyedge, nVert, nedge] =  midEDGEgen(nodeco, elnode, bdynde) ;
 
 GlobalV = zeros((nVert + nedge)*2, 1) ;
 
 bdynde=[];
 for i=1:size(nodeco,1)
     xx=nodeco(i,1);
     yy=nodeco(i,2);
     if xx*yy*(1-xx)*(1-yy)<1e-10
         bdynde = [bdynde;i];
     end
 end
 
% Get boundary info for ease in later boundary condition applications
% (Assumes if a node is dirichlet, then all components are dirichlet)
bdrysize = 2*size(bdynde,1);
bdrydof=[];
for ii=1:bdrysize/2
    bdrydof = [bdrydof; 2*bdynde(ii,1)-1; 2*bdynde(ii,1)];
end

if strcmp(pre_bas_type, 'CtsLin') == 1 
   GlobalP = zeros(nVert,2) ;

elseif strcmp(pre_bas_type, 'DCtsLin') == 1 
   ntri = size(elnode,1) ;
   GlobalP = zeros(3*ntri,2) ;
   GlobalP(:,2) = -1 ;
   mapN2GP = reshape(elnode(:,1:3).',[3*ntri 1]) ;
end
%%%

NVU = size(GlobalV,1) ;
NPU = size(GlobalP,1) ;
Sdim = NVU+NPU;
% create the matrices that will be involved in each timestep
%createConstantMatrices2;
%save constmats16bary MassMatrix StiffnessMatrix PressureMassMatrix PressureMatrix GradDivMatrix
load constmats16bary

% initialize u0
GlobalV = initVel(nodeco);
AllV = [(-dt+exp(-dt))*GlobalV,GlobalV];

gamma=1; 

errortable=[];
format shorte

% loop over timesteps
for t=dt:dt:T
    
    NSEsolve_BDF2
   
    AllV = [AllV(:,2),GlobalV];
    
    CalcErrDiv
        
    errortable=[errortable;VelL2Error,divL2,GradVelL2Error,PreL2Error]
    
    
end
    
%NSEsolveSV
%NSEsolveSV2
%NSEsolveSVYosida
    
