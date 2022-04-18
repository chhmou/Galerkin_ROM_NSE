function [localmat] = inner_prod_ten1_VDiv_ten1(triag_no, quad_rul, ...
   V_fun, ten1a_type, ten1b_type) 
%
% This function computes, for triangle triag_no, the integrals
% of the (vector) ten1a basis functions dot product with
% the (vector) Vfun Div ten1b basis functions.
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

% Evaluate the vector function V_fun at the quadrature points.
Vfun_vals = feval(V_fun, xy_pts, triag_no) ;

% Evaluate Basis Functions and their Gradients at quad. points.
[ten1a, Gradten1a] = feval(ten1a_type, quad_pts) ;
nbas1a = size(ten1a,1) ;

[ten1b, Gradten1b] = feval(ten1b_type, quad_pts) ;
nbas1b = size(ten1b,1) ;

% Do appropriate multiplies to get the true Gradients.
for iq = 1:nqpts
   Gradtrue1(:,:,iq) = Gradten1b(:,:,iq) * JInv ;
end


% Now to do the evaluations of the integrals.
for iq = 1:nqpts
   Gradtrue1(:,1,iq) = quad_wghts(iq) * Vfun_vals(1,iq) * Gradtrue1(:,1,iq) ;  
   Gradtrue1(:,2,iq) = quad_wghts(iq) * Vfun_vals(2,iq) * Gradtrue1(:,2,iq) ;  
end

VDivGrad(:,:) = Gradtrue1(:,1,:) + Gradtrue1(:,2,:) ;

mat1 = ten1a * VDivGrad.' ;

localmat = [ mat1, zeros(nbas1a, nbas1b)  ; ...
      zeros(nbas1a, nbas1b), mat1 ] ; 

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

