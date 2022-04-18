function [CQuadBVal, GradCQuadBVal] = CtsQuadBub(quad_pts)
%
% This function computes the values of the continuous
% quadratic basis functions, and of its gradient, at
% the quadrature points quad_pts --- on the reference triangle.
%

nqpt = size(quad_pts,1) ;

CQuadBVal(1,:) = (1.0 - quad_pts(:,1)' - quad_pts(:,2)').* ...
                (1.0 - 2*quad_pts(:,1)' - 2*quad_pts(:,2)') ;
CQuadBVal(2,:) = quad_pts(:,1)'.* (2*quad_pts(:,1)' - 1) ;
CQuadBVal(3,:) = quad_pts(:,2)'.* (2*quad_pts(:,2)' - 1) ;
CQuadBVal(4,:) = 4*quad_pts(:,1)'.* quad_pts(:,2)' ;
CQuadBVal(5,:) = 4*quad_pts(:,2)'.*( 1.0 - quad_pts(:,1)' - quad_pts(:,2)') ;
CQuadBVal(6,:) = 4*quad_pts(:,1)'.*( 1.0 - quad_pts(:,1)' - quad_pts(:,2)') ;
CQuadBVal(7,:) = 27*quad_pts(:,1)' .*quad_pts(:,2)' .*( 1.0 - quad_pts(:,1)' - quad_pts(:,2)') ;

GradCQuadBVal(1,1,:) = -3.0 + 4*quad_pts(:,1)' + 4*quad_pts(:,2)' ;
GradCQuadBVal(2,1,:) = 4*quad_pts(:,1)' - 1 ;
GradCQuadBVal(3,1,:) = zeros(1,nqpt) ;
GradCQuadBVal(4,1,:) = 4*quad_pts(:,2)' ;
GradCQuadBVal(5,1,:) = -4*quad_pts(:,2)' ;
GradCQuadBVal(6,1,:) = 4*(1.0 - 2*quad_pts(:,1)' - quad_pts(:,2)') ;
GradCQuadBVal(7,1,:) = 27*quad_pts(:,2)' .*( 1.0 - 2*quad_pts(:,1)' - quad_pts(:,2)') ;

GradCQuadBVal(1,2,:) = -3.0 + 4*quad_pts(:,1)' + 4*quad_pts(:,2)' ;
GradCQuadBVal(2,2,:) = zeros(1,nqpt) ;
GradCQuadBVal(3,2,:) = 4*quad_pts(:,2)' - 1 ;
GradCQuadBVal(4,2,:) = 4*quad_pts(:,1)' ;
GradCQuadBVal(5,2,:) = 4*(1.0 - quad_pts(:,1)' - 2*quad_pts(:,2)') ;
GradCQuadBVal(6,2,:) = -4*quad_pts(:,1)' ;
GradCQuadBVal(7,2,:) = 27*quad_pts(:,1)' .*( 1.0 - quad_pts(:,1)' - 2*quad_pts(:,2)') ;
