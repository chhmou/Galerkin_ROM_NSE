%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
%  This is the driver routine for the Navier-Stokes calculation
%
%

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  GlobalS  GlobalG
global dimTvel  dimTpre  dimTstr  dimTGrv 
global vel_bas_type  pre_bas_type  str_bas_type  Grv_bas_type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Define the problem and the solution method
%%

% We begin by specifying the type of approximation elements
% we are using for the velocity, pressure
% and the overall solution method.

%
vel_bas_type = 'CtsQuad' ;
pre_bas_type = 'CtsLin' ;


%sol_method = 'DefTen_SkSy_Nwt' ;


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Read in Grid.

R8x8 ;

%% For Driven Cavity -- Grid Generation
% xvec = [0 ; 0.5 ; 1] ;
% nxvec = [3 ; 3] ;
% yvec = [0 ; 0.5 ; 1] ;
% nyvec = [3 ; 3] ;
% [nodeco, elnode, eldegr, elnbgh, bdynde] =  gridrec(xvec, nxvec, yvec, nyvec) ;
% [bdynde] = bdysetupDC(nodeco, bdynde) ;
% [nodeco, elnode, bdynde, bdyedge, nVert, nedge] =  midEDGEgen(nodeco, elnode, bdynde) ;
%% End of Driven Cavity Grid Generation

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
   % Note: The velocity associated with node i is given by
   %       GlobalV(2*i-1,1) and GlobalV(2*i,1). Stored in
   %       GlobalV(:,2) is the boundary condition associates
   %       with node i.
   
   % Set up initial guess for the velocity and the boundary
   % condition variable. We assume that initVel imposes the
   % correct boundary conditions for the velocity and sets
   % the variable in GlobalV indicating if and what kind of  
   % bounday point it is.
   GlobalV(:,2) = -1 ;
   %
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
   
   % Set the boundary/initial data.
   tempV = initVel(nodeco) ;
   GlobalV(1:2:(nVert + nedge)*2,1) = tempV(1,:).' ;
   GlobalV(2:2:(nVert + nedge)*2,1) = tempV(2,:).' ;

     
elseif  strcmp(vel_bas_type, 'CtsQuadBub') == 1 
   ntri = size(elnode,1) ;
   GlobalV = zeros((nVert + nedge + ntri)*2, 2) ;
   % Note: The velocity associated with node i is given by
   %       GlobalV(2*i-1,1) and GlobalV(2*i,1). Stored in
   %       GlobalV(:,2) is the boundary condition associates
   %       with node i.
   
   % Set up initial guess for the velocity and the boundary
   % condition variable. We assume that initVel imposes the
   % correct boundary conditions for the velocity and sets
   % the variable in GlobalV indicating if and what kind of  
   % bounday point it is.
   GlobalV(:,2) = -1 ;
   %
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
   
   % Set the boundary/initial data. Note the coefficients of bubble fns are zero.
   tempV = initVel(nodeco) ;
   GlobalV(1:2:(nVert + nedge)*2,1) = tempV(1,:).' ;
   GlobalV(2:2:(nVert + nedge)*2,1) = tempV(2,:).' ;
   
elseif  strcmp(vel_bas_type, 'CtsLinBub') == 1 
   ntri = size(elnode,1) ;
   GlobalV = zeros((nVert + ntri)*2, 2) ;
   % Note: The velocity associated with node i is given by
   %       GlobalV(2*i-1,1) and GlobalV(2*i,1). Stored in
   %       GlobalV(:,2) is the boundary condition associates
   %       with node i.
   
   % Set up initial guess for the velocity and the boundary
   % condition variable. We assume that initVel imposes the
   % correct boundary conditions for the velocity and sets
   % the variable in GlobalV indicating if and what kind of  
   % bounday point it is.
   GlobalV(:,2) = -1 ;
   %
   for ibdy = 1 : size(bdynde,1)
      if ( bdynde(ibdy,2) ~= -1 )
         GlobalV( 2*bdynde(ibdy,1) - 1 , 2 ) = bdynde(ibdy,2) ;
         GlobalV( 2*bdynde(ibdy,1) , 2 ) = bdynde(ibdy,2) ;
      end
   end
   
   % Set the boundary/initial data. Note the coefficients of bubble fns are zero.
   tempV = initVel(nodeco(1:nVert,:)) ;
   GlobalV(1:2:2*nVert,1) = tempV(1,:).' ;
   GlobalV(2:2:2*nVert,1) = tempV(2,:).' ;

end


if strcmp(pre_bas_type, 'CtsLin') == 1 
   GlobalP = zeros(nVert,2) ;
   % Note: The pressure associated with traingle vertex i is
   %       given by GlobalP(i) 
   
   % Set up the initial guess for the Pressure and the boundary
   % condition variable. 
   % NOTE: Through there is generally not a boundary condition for
   %       the pressure we include this here to enable the pressure
   %       to be specified at a point to provide uniqueness for the
   %       pressure.
   % NOTE NOTE: As we normalize the pressure to zero at this point
   %            this is already taken care of by the definition of GlobalP
   GlobalP(:,2) = -1 ;
   %
   for ibdy = 1 : size(bdynde,1)
      if ( bdynde(ibdy,3) == 0 )
         GlobalP( bdynde(ibdy,1) , 2 ) = 0 ;
      end
   end
   
elseif strcmp(pre_bas_type, 'DCtsLin') == 1 
   ntri = size(elnode,1) ;
   GlobalP = zeros(3*ntri,2) ;
   % Note: The pressure associated with triangle vertex i is
   %       given by GlobalP(i) 
   
   % Set up the initial guess for the Pressure and the boundary
   % condition variable. 
   % NOTE: Through there is generally not a boundary condition for
   %       the pressure we include this here to enable the pressure
   %       to be specified at a point to provide uniqueness for the
   %       pressure.
   % NOTE NOTE: As we normalize the pressure to zero at this point
   %            this is already taken care of by the definition of GlobalP
   GlobalP(:,2) = -1 ;
   %
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


%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Now to solve the problem !  

%if strcmp(sol_method, 'DefTen_SkSy_Nwt') == 1
%   simplenwt
%   deftennwt
   deftensksynwt
%end


%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Postprocess the approximation.

% Call a subroutine to graph the various quantities.


      
