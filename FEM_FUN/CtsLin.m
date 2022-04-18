function [CLinVal, GradCLinVal] = CtsLin(quad_pts)
%
% This function computes the values of the continuous
% linear basis functions, and of its gradient, at
% the quadrature points quad_pts --- on the reference triangle.
%

nqpt = size(quad_pts,1) ;

CLinVal(1,:) = 1.0 - quad_pts(:,1)' - quad_pts(:,2)' ;
CLinVal(2,:) = quad_pts(:,1)' ;
CLinVal(3,:) = quad_pts(:,2)' ;


GradCLinVal(1,1,:) = -1*ones(1,nqpt) ;
GradCLinVal(2,1,:) =  1*ones(1,nqpt) ;
GradCLinVal(3,1,:) =  zeros(1,nqpt) ;

GradCLinVal(1,2,:) = -1*ones(1,nqpt) ;
GradCLinVal(2,2,:) = zeros(1,nqpt) ;
GradCLinVal(3,2,:) = 1*ones(1,nqpt) ;
