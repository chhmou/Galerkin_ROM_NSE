%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDF3 + algebraic splitting


clear
format long

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge GlobalVnew GlobalPSV
global GlobalV  GlobalP  nu vel_bas_type  pre_bas_type problem sig sumV





vel_bas_type = 'CtsQuad' ;  % no choice at the moment
pre_bas_type = 'DCtsLin' ;  % could be DCtsLin

nu=1/20;
t=0;
gammas=[1,10,100,1000,10000];
sig=1;



    
problem = 'Step'

[nodeco, elnode, bdynde] = FFEM2ErvGrd2('SVgridStep.msh') ;


nbdy = size(bdynde,1) ;
bdynde = [bdynde , zeros(nbdy , 1), -1*ones(nbdy , 2) ] ;
[nodeco, elnode, bdynde, bdyedge, nVert, nedge] =  midEDGEgen(nodeco, elnode, bdynde) ;
 
 GlobalV = zeros((nVert + nedge)*2, 1) ;
 
 bdynde=[];
 for i=1:size(nodeco,1)
     xx=nodeco(i,1);
     yy=nodeco(i,2);
     if xx*yy*(40-xx)*(10-yy)<1e-10 || ( xx>4.9999 && xx<6.0001 && yy<1.00001 )
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
%save stepmatsbary MassMatrix StiffnessMatrix PressureMassMatrix PressureMatrix GradDivMatrix
load stepmatsbary

% initialize for U in Udotgradu and boundary conditions
GlobalV = initVel(nodeco);

% solve the coupled problem
gamma=100;
NSEsolveSV
% store the solution
GlobalVSV = GlobalV;
GlobalPSV = GlobalP;
%plotStep2;

%   gamma=10000;
%   GlobalV = initVel(nodeco);
%   NSEsolveSV2
% 

errortable=[];

for gindex=1:5

    gamma=gammas(gindex);



    % % initialize for U in Udotgradu and boundary conditions
    % GlobalV = initVel(nodeco);
    % NSEsolveSV2

    GlobalV = initVel(nodeco);
    NSEsolveSVYosida
    %plotStep2

    CalcErrDiv
    GlobalV = GlobalV - GlobalVSV;
    % find || grad GlobalV || = || u_h ||
    CalcErrDiv2

    format shorte
    if gindex<3
        row = [H1norm, 0, divL2, 0];
    else
        row = [H1norm, log(errortable(gindex-1,1)/H1norm)/log(10),...
            divL2, log(errortable(gindex-1,3)/divL2)/log(10)];
    end
    errortable=[errortable;row]
     [preErr1, preErr2, preErr3, preErr4]
end
% 
% store the solution
GlobalV = GlobalVSV;
CalcErrDiv
errortable = [errortable; 0,0,divL2,0] 
%       
