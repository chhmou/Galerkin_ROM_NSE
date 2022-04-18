function [localmat] = Identity(xyz_pts, tet_no, Rxyz_pts, Jinv) 
%
% This function computes, the current approximation for the velocity
% at the requested xy_pts points in triangle triag_no.
%  The matrix of values is returned in localmat.
%  
%

%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global nodeco  elnode  bdynde  bdyedge  nVert  nedge
global GlobalV  GlobalP  GlobalS  GlobalG
global dimTvel  dimTpre  dimTstr  dimTGrv 
global vel_bas_type  pre_bas_type  str_bas_type  Grv_bas_type


localmat(1,1,:) = 1 ;
localmat(1,2,:) = 0 ;
localmat(2,1,:) = 0 ;
localmat(2,2,:) = 1 ;
