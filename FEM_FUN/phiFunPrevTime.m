function [localmat] = phiFunPrevTime(xy_pts, triag_no) 
%
% This function computes, the current approximation for the velocity
% at the requested xy_pts points in triangle triag_no.
%  The matrix of values is returned in localmat.
%  
%

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  PastV PastPhi
global vel_bas_type  pre_bas_type


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

   Vstart = [2*(elnode(triag_no,1:3) - 1) + 1 , 2*(elnode(triag_no,4:6) + nVert - 1) + 1 ] ;
   Vel1 = PastPhi([Vstart], end) ;
   Vstart = Vstart + 1 ;
   Vel2 = PastPhi([Vstart], end) ;
%Vstart = [2*(elnode(triag_no,1:6) - 1) + 1] ;
%Vel1 = PastV([Vstart], end) ;
%Vstart = Vstart + 1 ;
%Vel2 = PastV([Vstart], end) ;
      
% Multiply the coefficients by the basis values.
localmat(1,:) = Vel1.' * basvals ;
localmat(2,:) = Vel2.' * basvals ;
