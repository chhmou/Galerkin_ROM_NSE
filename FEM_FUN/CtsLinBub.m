function [CLinBVal, GradCLinBVal] = CtsLinBub(quad_pts)
%
% This function computes the values of the continuous
% linear basis functions, and of its gradient, at
% the quadrature points quad_pts --- on the reference triangle.
%

nqpt = size(quad_pts,1) ;

CLinBVal(1,:) = 1.0 - quad_pts(:,1)' - quad_pts(:,2)' ;
CLinBVal(2,:) = quad_pts(:,1)' ;
CLinBVal(3,:) = quad_pts(:,2)' ;
CLinBVal(4,:) = 27*quad_pts(:,1)' .*quad_pts(:,2)' .*( 1.0 - quad_pts(:,1)' - quad_pts(:,2)') ;

GradCLinBVal(1,1,:) = -1*ones(1,nqpt) ;
GradCLinBVal(2,1,:) =  1*ones(1,nqpt) ;
GradCLinBVal(3,1,:) =  zeros(1,nqpt) ;
GradCLinBVal(4,1,:) = 27*quad_pts(:,2)' .*( 1.0 - 2*quad_pts(:,1)' - quad_pts(:,2)') ;

GradCLinBVal(1,2,:) = -1*ones(1,nqpt) ;
GradCLinBVal(2,2,:) = zeros(1,nqpt) ;
GradCLinBVal(3,2,:) = 1*ones(1,nqpt) ;
GradCLinBVal(4,2,:) = 27*quad_pts(:,1)' .*( 1.0 - quad_pts(:,1)' - 2*quad_pts(:,2)') ;
