function [velocity2] = initVel2(xy_pts,time) 


xpts = xy_pts(:,1).' ;
ypts = xy_pts(:,2).' ;
npts = size(xy_pts,1) ;

velocity = zeros(2,npts) ;



    velocity(1,1:npts) = sin(pi*xpts).*sin(pi*ypts)*exp(time);
    velocity(2,1:npts) = cos(pi*xpts).*cos(pi*ypts)*exp(time);

velocity2 = zeros(2*npts,1);
velocity2(1:2:end) = velocity(1,:)' ;
velocity2(2:2:end) = velocity(2,:)' ;
