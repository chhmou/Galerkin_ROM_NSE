function [localmat] = extrapvelGradvelPrev(xy_pts, triag_no) 
%
% This function computes, the current approximation for the 
% product of velocity dot gradient velocity at the requested 
% xy_pts points in triangle triag_no.
%  The matrix of values is returned in localmat.
%  
%

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  AllV DNbarOfUNEW PastV
global dimTvel  dimTpre  dimTstr  dimTGrv 
global vel_bas_type  pre_bas_type  str_bas_type  Grv_bas_type


%%%%
npts = size(xy_pts,1) ;

% Get the velocity and Gradient of velocity values
Vfun_vals = extrapvelFun(xy_pts, triag_no);
Grad_vals = GradvelFunPrevTime(xy_pts, triag_no) ;

%% In MATLAB multiplication of other than matrices is not defined so we 
%% need to make some provisions.

avec1x(1:npts) = Grad_vals(1,1,1:npts) ;
avec1y(1:npts) = Grad_vals(2,1,1:npts) ;
avec2x(1:npts) = Grad_vals(1,2,1:npts) ;
avec2y(1:npts) = Grad_vals(2,2,1:npts) ;

   
localmat(1,:) = Vfun_vals(1,:) .* avec1x + Vfun_vals(2,:) .* avec1y ;
localmat(2,:) = Vfun_vals(1,:) .* avec2x + Vfun_vals(2,:) .* avec2y ;






