function [localmat] = etafun(xy_pts, triag_no) 
%
% This function computes, the values for etafun -- the
% viscosity term which appears in the momentum equation
% at the requested xy_pts points in triangle triag_no.
%  The vector of values is returned in localmat.
%  
%

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  GlobalS  GlobalG nu
global dimTvel  dimTpre  dimTstr  dimTGrv 
global vel_bas_type  pre_bas_type  str_bas_type  Grv_bas_type


%%%%
npts = size(xy_pts,1) ;

%% localmat is a vector of values
localmat = nu * ones(1,npts) ;





