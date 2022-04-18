function [localmat] = ptrue(xy_pts, triag_no,time) 
%
% This function defines the true value for the pressure.
% If the exact value for the pressure is not know, entering
% zero will result in the norm of the pressure being calculated.
%

global t

%%%%
xpts = xy_pts(:,1).' ;
ypts = xy_pts(:,2).' ;
npts = size(xy_pts,1) ;

localmat = zeros(1,npts) ;

%localmat = xpts+ypts ;

%localmat = 100*sin(2*(xpts+ypts));

%localmat(1:npts) = 2*xpts + 3*ypts - 2 ;

%localmat(1:npts) = ypts.^2 .* xpts.^2 + (xpts - ypts - 1).* sin(pi*xpts/2) ;

% True pressue for the test example for the Scott-Vogelius element
%localmat(1:npts) = 2.0* cos(xpts) .*sin(ypts) - 2.0* sin(1.0) *(1.0 - cos(1.0)) ;

%localmat(1:npts) = xpts + ypts + 1/2*(1+.01*thetime)^2*( (cos(ypts)).^2 + (sin(xpts)).^2 ) - 1 - 1/2*(1+.01*thetime)^2*(cos(1))^2;


%localmat(1:npts) = cos(xpts).*sin(ypts+t);

localmat(1:npts) = sin(xpts+ypts)*exp(2*t);

