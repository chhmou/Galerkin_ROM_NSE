function [localmat] = inner_prod_ten1_Vec(triag_no, quad_rul, ...
   Vec_fun, ten1_type) 
%
% This function computes, for triangle triag_no, the integrals
% of the (vector) ten1a basis functions times 
% the vector function Vec_fun. 
% The vector of values is returned in localmat.
%  
%


%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  GlobalS  GlobalG
global dimTvel  dimTpre  dimTstr  dimTGrv 
global vel_bas_type  pre_bas_type  str_bas_type  Grv_bas_type


% Description of triangle.
cotri(1:3,1) = nodeco(elnode(triag_no, 1:3), 1) ;
cotri(1:3,2) = nodeco(elnode(triag_no, 1:3), 2) ;
    
Jmat = [(cotri(2,1) - cotri(1,1)), (cotri(3,1) - cotri(1,1)) ; ...
        (cotri(2,2) - cotri(1,2)) , (cotri(3,2) - cotri(1,2)) ] ;
detJ = abs(Jmat(1,1)*Jmat(2,2) - Jmat(1,2)*Jmat(2,1));
JInv = inv(Jmat) ;


% Evaluation of quadrature points and quadrature weights.
[quad_pts, quad_wghts] = feval(quad_rul) ;
nqpts = size(quad_pts,1) ;

% Adjust points and weights to account for size of true triangle.
xy_pts = ( Jmat * quad_pts.' ).' ;
xy_pts(:,1) = cotri(1,1) + xy_pts(:,1) ;
xy_pts(:,2) = cotri(1,2) + xy_pts(:,2) ;
quad_wghts = detJ * quad_wghts ;

% Evaluate the Vector function multiplier at the quadrature points.
Vfun_vals = feval(Vec_fun, xy_pts, triag_no) ;

% Evaluate Basis Functions and their Gradients at quad. points.
[ten1, Gradten1] = feval(ten1_type, quad_pts) ;
nbas1 = size(ten1,1) ;


% Now to do the evaluations of the integrals.
for iq = 1:nqpts
   ten1(:,iq) = quad_wghts(iq) * ten1(:,iq) ;  
end

mat1 = ten1 * Vfun_vals(1,:).' ;
mat2 = ten1 * Vfun_vals(2,:).' ;

% Note: localmat is a vector
localmat = [ mat1 ; mat2 ] ; 

%%
%% We now must reorder the matrix so that the numbering of the unknowns
%% is consistent with that of sequentially numbering the unknowns at each
%% node ---- as opposed to sequentially numbering the unknowns relative to
%% the basis representing the unknown.

[rdim cdim] = size(localmat) ;
rowperm = [ ] ; colperm = [ ] ; 
for ir = 1:nbas1
   rowperm = [rowperm ir:nbas1:rdim] ;
end
localmat = localmat(rowperm) ;

