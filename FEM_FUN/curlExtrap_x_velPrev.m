function [localmat] = curlExtrap_x_velPrev(xy_pts, triag_no) 
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
curlExtrap_as_scalar = curlExtrap_2d(xy_pts, triag_no);
velPrevFunVals = velFunPrevTime(xy_pts, triag_no) ;
   
localmat(1,:) = - curlExtrap_as_scalar .* velPrevFunVals(2,:);
localmat(2,:) =   curlExtrap_as_scalar .* velPrevFunVals(1,:);





