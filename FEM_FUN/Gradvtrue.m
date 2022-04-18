function [localmat] = Gradvtrue(xy_pts, triag_no,time) 
%
% This function defines the true value for the velocity.
% If the exact value for the velocity is not know, entering
% zero will result in the norm of the velocity being calculated.
%

global nu problem t

%%%%
xpts = xy_pts(:,1).' ;
ypts = xy_pts(:,2).' ;
npts = size(xy_pts,1) ;

localmat = zeros(2,2,npts) ;

 if strcmp(problem,'SplitTest1')

 localmat(1,1,1:npts) = pi*cos(pi*xpts).*sin(pi*ypts);
 localmat(2,1,1:npts) = pi*sin(pi*xpts).*cos(pi*ypts);
 localmat(1,2,1:npts) = -pi*sin(pi*xpts).*cos(pi*ypts); 
 localmat(2,2,1:npts) = -pi*cos(pi*xpts).*sin(pi*ypts);
 localmat = exp(t)*localmat;
 elseif strcmp(problem,'YosidaTest')
    localmat(1,1,1:npts) = 0;
    localmat(2,1,1:npts) = -sin(ypts);
    localmat(1,2,1:npts) = cos(xpts); 
    localmat(2,2,1:npts) = 0; 
 elseif strcmp(problem,'YosidaTest2')
    localmat(1,1,1:npts) = 0;
    localmat(2,1,1:npts) = -sin(ypts);
    localmat(1,2,1:npts) = cos(xpts); 
    localmat(2,2,1:npts) = 0;
    localmat = localmat * (t+exp(t));
 elseif strcmp(problem,'Step')
     
 else
     error('problem invalid in gradvtrue')
 end