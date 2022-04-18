function [localmat] = inner_prod_ten0_Div_ten1(triag_no, quad_rul, ...
   scal_fun, ten0_type, ten1_type) 
%
% This function computes, for triangle triag_no, the integrals
% of the (scalar) ten0 basis functions times the divergence of
% the (vector) ten1 basis functions multiplied by the scalar
% function scal_fun. The matrix of values is returned in localmat.
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
[ten0, Gradten0] = feval(ten0_type, quad_pts) ;
nbas0 = size(ten0,1) ;

[ten1, Gradten1] = feval(ten1_type, quad_pts) ;
nbas1 = size(ten1,1) ;

% Do appropriate multiplies to get the true Gradients.
for iq = 1:nqpts
   Gradtrue(:,:,iq) = Gradten1(:,:,iq) * JInv ;
end


% Now to do the evaluations of the integrals.
for iq = 1:nqpts
%   DivGradtrue(:,iq) = Gradtrue(:,:,iq) * [1 ; 1] ;
   ten0(:,iq) = quad_wghts(iq) * sfun_vals(iq) * ten0(:,iq) ;  
end

mat1(:,:) = Gradtrue(:,1,:) ;
mat2(:,:) = Gradtrue(:,2,:) ;
localmat = [ten0 * mat1.' , ten0 * mat2.'] ;

%%
%% We now must reorder the matrix so that the numbering of the unknowns
%% is consistent with that of sequentially numbering the unknowns at each
%% node ---- as opposed to sequentially numbering the unknowns relative to
%% the basis representing the unknown.

[rdim cdim] = size(localmat) ;
rowperm = [ ] ; colperm = [ ] ; 
for ic = 1:nbas1
   colperm = [colperm ic:nbas1:cdim] ;
end
localmat = localmat(: , colperm) ;


