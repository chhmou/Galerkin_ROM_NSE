function [localmat] = curlVelocity_2d(xy_pts, triag_no) 
%
% note this returns a scalar!!!  WE are in 2d, so the curl is actually a 3d
% vector but only nonzero component is 3rd component.  If we returned the
% whole 3d vectro, this wouldn't fit well with code.

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  PastV ExtrapTerm
global vel_bas_type  pre_bas_type  str_bas_type  Grv_bas_type


%%%%
% We firstly need to determine the corresponding points on the
% reference triangle.

% Description of triangle.
cotri(1:3,1) = nodeco(elnode(triag_no, 1:3), 1) ;
cotri(1:3,2) = nodeco(elnode(triag_no, 1:3), 2) ;
    
Jmat = [(cotri(2,1) - cotri(1,1)), (cotri(3,1) - cotri(1,1)) ; ...
        (cotri(2,2) - cotri(1,2)) , (cotri(3,2) - cotri(1,2)) ] ;
detJ = abs(Jmat(1,1)*Jmat(2,2) - Jmat(1,2)*Jmat(2,1));
JInv = inv(Jmat) ;


% Map the points to the reference triangle.
Rxy_pts(:,1) = xy_pts(:,1) - cotri(1,1) ; 
Rxy_pts(:,2) = xy_pts(:,2) - cotri(1,2) ;
Rxy_pts = Rxy_pts * JInv.' ;

% Evaluate Basis Functions and their Gradients at requested points.
[basvals, Gradten1] = feval(vel_bas_type, Rxy_pts) ;
nqpts = size(basvals,2) ;

% Do appropriate multiplies to get the true Gradients.
for iq = 1:nqpts
   Gradtrue(:,:,iq) = Gradten1(:,:,iq) * JInv ;
end

% Extract from the Global solution vector the coefficient values for
% the basis functions.
   Vstart = [2*(elnode(triag_no,1:3) - 1) + 1 , 2*(elnode(triag_no,4:6) + nVert - 1) + 1 ] ;
   Vel1 = GlobalV([Vstart],1) ;
   Vstart = Vstart + 1 ;
   Vel2 = GlobalV([Vstart],1) ;

%% In MATLAB multiplication of other than matrices is not defined so we 
%% need to make some provisions.

Gradx(:,:) = Gradtrue(:,1,:) ;
Grady(:,:) = Gradtrue(:,2,:) ;
   
%% Multiply the coefficients by the basis values.
localmat = Vel2.' * Gradx - Vel1.' * Grady;







