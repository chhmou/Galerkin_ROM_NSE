function [localmat] = inner_prod_Grad_ten1T_Grad_ten1(triag_no, quad_rul, ...
   scal_fun, ten1a_type, ten1b_type) 
%
% This function computes, for triangle triag_no, the integrals
% of the TRANSPOSE of the gradient of (vector) ten1a basis functions times 
% the gradient of (vector) ten1b basis functions 
% multiplied by the scalar function scal_fun. 
%
% **** Identical to inner_prod_Grad_ten1_Grad_ten1T  ******
% 
% The matrix of values is returned in localmat.
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

% Evaluate the scalar multiplier at the quadrature points.
sfun_vals = feval(scal_fun, xy_pts, triag_no) ;

% Evaluate Basis Functions and their Gradients at quad. points.
[ten1a, Gradten1a] = feval(ten1a_type, quad_pts) ;
nbas1a = size(ten1a,1) ;

[ten1b, Gradten1b] = feval(ten1b_type, quad_pts) ;
nbas1b = size(ten1b,1) ;

% Do appropriate multiplies to get the true Gradients.
for iq = 1:nqpts
   Gradtrue1(:,:,iq) = Gradten1a(:,:,iq) * JInv ;
   Gradtrue2(:,:,iq) = Gradten1b(:,:,iq) * JInv ;
end


% Now to do the evaluations of the integrals.
for iq = 1:nqpts
   Gradtrue1(:,:,iq) = quad_wghts(iq) * sfun_vals(iq) * Gradtrue1(:,:,iq) ;  
end

tempM11(:,:) = Gradtrue1(:,1,:) ;
tempM12(:,:) = Gradtrue1(:,2,:) ;
tempM21(:,:) = Gradtrue2(:,1,:) ;
tempM22(:,:) = Gradtrue2(:,2,:) ;
mat11 = tempM11 * tempM21.' ;
mat12 = tempM12 * tempM21.' ;
mat21 = tempM11 * tempM22.' ;
mat22 = tempM12 * tempM22.' ;

localmat = [ mat11 , mat12 ; ...
             mat21 , mat22 ] ;

%%
%% We now must reorder the matrix so that the numbering of the unknowns
%% is consistent with that of sequentially numbering the unknowns at each
%% node ---- as opposed to sequentially numbering the unknowns relative to
%% the basis representing the unknown.

[rdim cdim] = size(localmat) ;
rowperm = [ ] ; colperm = [ ] ; 
for ir = 1:nbas1a
   rowperm = [rowperm ir:nbas1a:rdim] ;
end
for ic = 1:nbas1b
   colperm = [colperm ic:nbas1b:cdim] ;
end
localmat = localmat(rowperm , colperm) ;



