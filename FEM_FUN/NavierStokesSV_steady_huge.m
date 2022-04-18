% Testing whether (P2,0) with barycenter mesh and grad-div is any good

% 8/6/10
% Leo Rebholz

clear

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  GlobalS  GlobalG nu
global dimTvel  dimTpre  dimTstr  dimTGrv 
global vel_bas_type  pre_bas_type  str_bas_type  Grv_bas_type


%
vel_bas_type = 'CtsQuad' ;
pre_bas_type = 'DCtsLin' ;

nu=1/100;
gd_parameter=1;

% Enter the grid file
GenSVgrd('SQ8.msh');
[nodeco, elnode, bdynde] = FFEM2ErvGrd2('SVgrid.msh') ;

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
createConstantMatrices2;

%% Now go through the time steps and solve!
    

%% Solve w/ (P2,0) elements
NSEsolve_nonlinear_steady_SV_SS_huge


% 
% 
% % reinitialize problem and run rotational form
% load SETUP
% NSEsolve_nonlinear_steady_SV_rot
% rot_vel = GlobalV(:,1);
% rot_pre = GlobalP(:,1);
% CalcErrDiv
% 
% 
% % reinitalize problem and run skew-symmetric form
% load SETUP
% NSEsolve_nonlinear_steady_SV_SS
% ss_vel = GlobalV(:,1);
% ss_pre = GlobalP(:,1);
% CalcErrDiv
% 




              
