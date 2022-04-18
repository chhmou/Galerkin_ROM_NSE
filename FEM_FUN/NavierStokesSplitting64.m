%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDF3 + algebraic splitting


clear
clc

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge t PastV
global GlobalV  GlobalP  nu vel_bas_type  pre_bas_type problem




% 
% vel_bas_type = 'CtsQuad' ;  % no choice at the moment
% pre_bas_type = 'DCtsLin' ;  % could be DCtsLin
% 
% nu=1.0;
% T=4;
% dt=T/4;
% t=0;
% numTimeSteps=ceil(T/dt);
% 
% gamma=0;
% 
% problem = 'SplitTest1'
% 
% meshname = 'SquareMeshU64.msh';
% 
% 
% %%%
% GenSVgrd(meshname);
% [nodeco, elnode, bdynde] = FFEM2ErvGrd2('SVgrid.msh') ;
% 
% nbdy = size(bdynde,1) ;
% bdynde = [bdynde , zeros(nbdy , 1), -1*ones(nbdy , 2) ] ;
% [nodeco, elnode, bdynde, bdyedge, nVert, nedge] =  midEDGEgen(nodeco, elnode, bdynde) ;
%  
%  GlobalV = zeros((nVert + nedge)*2, 1) ;
%  
%  bdynde=[];
%  for i=1:size(nodeco,1)
%      xx=nodeco(i,1);
%      yy=nodeco(i,2);
%      if xx*yy*(1-xx)*(1-yy)<1e-10
%          bdynde = [bdynde;i];
%      end
%  end
%  
% % Get boundary info for ease in later boundary condition applications
% % (Assumes if a node is dirichlet, then all components are dirichlet)
% bdrysize = 2*size(bdynde,1);
% bdrydof=[];
% for ii=1:bdrysize/2
%     bdrydof = [bdrydof; 2*bdynde(ii,1)-1; 2*bdynde(ii,1)];
% end
% 
% if strcmp(pre_bas_type, 'CtsLin') == 1 
%    GlobalP = zeros(nVert,2) ;
% 
% elseif strcmp(pre_bas_type, 'DCtsLin') == 1 
%    ntri = size(elnode,1) ;
%    GlobalP = zeros(3*ntri,2) ;
%    GlobalP(:,2) = -1 ;
%    mapN2GP = reshape(elnode(:,1:3).',[3*ntri 1]) ;
% end
% %%%
% 

load h64P2P1dcData
dt=T/16;
t=0;
numTimeSteps=ceil(T/dt);


% % We are using BDF3, so I need 3 initial timesteps - take them to be true
% % solution
 PastV = [initVel2(nodeco,0),initVel2(nodeco,dt),initVel2(nodeco,2*dt)];
 GlobalV = PastV(:,end);
 PastP = [GlobalP(:,1),GlobalP(:,1),GlobalP(:,1)];

% % create the matrices that will be involved in each timestep
% createConstantMatrices2;
 

% t=0;
% GlobalV = PastV(:,end-2);
% t=dt;
% GlobalV = PastV(:,end-1);
% CalcErrDiv
% 
% t=2*dt;
% GlobalV = PastV(:,end);
% CalcErrDiv
% 
% save h64P2P1dcData


% Now go through the time steps and solve!
for jj=3:numTimeSteps

    t = jj*dt;
    display('**********************')
    NSEsolve_BDF3AB3
            
    % store the solution
    PastV = [PastV(:,2:3),GlobalV];
    PastP = [PastP(:,2:3),GlobalP];
    %CalcErrDiv
   
end
CalcErrDiv   
NVU

      
