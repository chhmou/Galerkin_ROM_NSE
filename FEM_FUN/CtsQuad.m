function [CQuadVal, GradCQuadVal] = CtsQuad(quad_pts)
%
% This function computes the values of the continuous
% quadratic basis functions, and of its gradient, at
% the quadrature points quad_pts --- on the reference triangle.
%

nqpt = size(quad_pts,1) ;

CQuadVal(1,:) = (1.0 - quad_pts(:,1)' - quad_pts(:,2)').* ...
                (1.0 - 2*quad_pts(:,1)' - 2*quad_pts(:,2)') ;
CQuadVal(2,:) = quad_pts(:,1)'.* (2*quad_pts(:,1)' - 1) ;
CQuadVal(3,:) = quad_pts(:,2)'.* (2*quad_pts(:,2)' - 1) ;
CQuadVal(4,:) = 4*quad_pts(:,1)'.* quad_pts(:,2)' ;
CQuadVal(5,:) = 4*quad_pts(:,2)'.*( 1.0 - quad_pts(:,1)' - quad_pts(:,2)') ;
CQuadVal(6,:) = 4*quad_pts(:,1)'.*( 1.0 - quad_pts(:,1)' - quad_pts(:,2)') ;

GradCQuadVal(1,1,:) = -3.0 + 4*quad_pts(:,1)' + 4*quad_pts(:,2)' ;
GradCQuadVal(2,1,:) = 4*quad_pts(:,1)' - 1 ;
GradCQuadVal(3,1,:) = zeros(1,nqpt) ;
GradCQuadVal(4,1,:) = 4*quad_pts(:,2)' ;
GradCQuadVal(5,1,:) = -4*quad_pts(:,2)' ;
GradCQuadVal(6,1,:) = 4*(1.0 - 2*quad_pts(:,1)' - quad_pts(:,2)') ;

GradCQuadVal(1,2,:) = -3.0 + 4*quad_pts(:,1)' + 4*quad_pts(:,2)' ;
GradCQuadVal(2,2,:) = zeros(1,nqpt) ;
GradCQuadVal(3,2,:) = 4*quad_pts(:,2)' - 1 ;
GradCQuadVal(4,2,:) = 4*quad_pts(:,1)' ;
GradCQuadVal(5,2,:) = 4*(1.0 - quad_pts(:,1)' - 2*quad_pts(:,2)') ;
GradCQuadVal(6,2,:) = -4*quad_pts(:,1)' ;
