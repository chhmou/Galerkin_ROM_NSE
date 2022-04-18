%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear
format long

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge GlobalVnew
global GlobalV  GlobalP  nu vel_bas_type  pre_bas_type problem sig
global GlobalPSV sumV t AllV dt




vel_bas_type = 'CtsQuad' ;  
pre_bas_type = 'CtsLin' ;  

 problem = 'Step'
 nu=1/600;
 t=0;
 T=40;
 dt=0.2;
 meshname = 'Stepmesh1.msh';
% [nodeco, elnode, bdynde] = FFEM2ErvGrd2('Stepmesh1.msh') ;
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
     if strcmp('YosidaTest2',problem)
         if xx*yy*(1-xx)*(1-yy)<1e-10
             bdynde = [bdynde;i];
         end
     else
         % step
         if abs( xx*yy*(10-yy)*(40-xx) )<1e-10 || (xx>4.9999 && xx<6.00001 && yy<1.00001)
             bdynde = [bdynde;i];
         end
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
createConstantMatrices2;
%save constmats16bary MassMatrix StiffnessMatrix PressureMassMatrix PressureMatrix GradDivMatrix

% for YosidaTest2 with bary16 mesh and SV elts
%load constmats16bary

% initialize u0
GlobalV = initVel(nodeco);
AllV = [GlobalV,GlobalV];

gamma=1.0; 

errortable=[];
format shorte

% first timestep, use usual BDF2 w/out approximation
t=dt;
NSEsolve_BDF2_noYosida
AllV = [AllV(:,2),GlobalV];
%CalcErrDiv
%errortable=[errortable;VelL2Error,divL2,GradVelL2Error,PreL2Error]

index=0;
% loop over timesteps 2 to M, and solve for updates
for t=2*dt:dt:T
    index=index+1;
     display(['t=' num2str(t) ', BDF2 Yosida, G=1 TH, barymesh, dt=0.2' ])
     % usual BDF2
%     NSEsolve_BDF2_noYosida
    
     % BDF2, with Yosida
     NSEsolve_BDF2_Yosida
    
    % BDF2, with Yosida + pressure correction
%    NSEsolve_BDF2_YosidaCorrected
  
       % BDF2 solving for updates, Yosida
%      NSEsolve_BDF2updates_Yosida
      
   
       % BDF2 solving for updates, Yosida + pressure correction
%      NSEsolve_BDF2updates_YosidaCorrected
 
% if mod(index,10)==0
%     NSEsolve_BDF2_YosidaCorrected
% else
%      NSEsolve_BDF2updates_Yosida
% end

    AllV = [AllV(:,2),GlobalV];
    
%    CalcErrDiv
%    divL2
        
end




