%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDF3 + algebraic splitting

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  nu vel_bas_type  pre_bas_type problem

vel_bas_type = 'CtsQuad' ;  % no choice at the moment
pre_bas_type = 'CtsLin' ;  % could be DCtsLin

nu=1.0;
T=1;
dt=T/4;
numTimeSteps=ceil(T/dt);

problem = 'SplitTest1'

meshname = 'SQ16.msh';

createMeshandElementStructure;
% nodeco is list of coordinates of nodes
% elnode is list of nodes for each element
% GlobalV and GlobalP are initialized to be zero
% bdynde is the list of boundary nodes, and bdrydof is their global dof #'s

% We are using BDF3, so I need 3 initial timesteps - take them to be true
% solution
PastV = [initVel2(nodeco,0),initVel2(nodeco,dt),initVel2(nodeco,2*dt)];
PastP = [GlobalP(:,1),GlobalP(:,1),GlobalP(:,1)];

% create the matrices that will be involved in each timestep
CreateConstantMatrices2;


% Now go through the time steps and solve!
for jj=3:numTimeSteps

    t = jj*dt;
    
    NSEsolve_nonlinear_SV
            
    % store the solution
    PastV = [PastV(:,2:3),GlobalV(:,1)];
    PastP = [PastP(:,2:3),GlobalP(:,1)];

   
end

CalcErrDiv

      
