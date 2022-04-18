function [velocity] = initVel(xy_pts,bdryelts) 

%% This function sets up the boundary condtion for the velocity
%% and gives an initial approximation to the velocity throughout
%% the domain.
%% For the channel flow problem we assume a parabolic velocity
%% inflow profile with rate of flow equal to 1.
%% As the initial guess throughout the domain we assume a parabolic
%% profile and adjust the "amplitude" so that the conservation of
%% mass is preserved.

global problem

xpts = xy_pts(:,1).' ;
ypts = xy_pts(:,2).' ;
npts = size(xy_pts,1) ;

velocity = zeros(2,npts) ;

if strcmp(problem,'graddivtest')

    
% We have to handle the problem of \int_{edge} u \cdot n = \int_{edge} uh \cdot n
% Loop over the boundary edges
%for ii=1:size(bdryelts,1)
    elnumber = bdryelts(ii,1);
    activeedge = bdryelts(ii,2);
    
    
    
%end


    velocity(1,1:npts) = xpts.*sin(ypts);
    velocity(2,1:npts) = cos(ypts);
elseif strcmp(problem,'graddivtest2')
    velocity(1,1:npts) = sin(ypts);
    velocity(2,1:npts) = cos(xpts);
elseif strcmp(problem,'YosidaTest') || strcmp(problem,'YosidaTest2')
    velocity(1,1:npts) = cos(ypts);
    velocity(2,1:npts) = sin(xpts);
    
    velocity2 = zeros(2*npts,1);
    velocity2(1:2:end) = velocity(1,:)' ;
    velocity2(2:2:end) = velocity(2,:)' ;
    velocity=velocity2;
elseif strcmp(problem,'Step')
    velocity(1,1:npts) = ypts.*(10-ypts)/25;
    velocity(2,1:npts) = 0;
    
    for i=1:npts
        if xpts(i)>4.999 && xpts(i)<6.0001  && ypts(i)<1.0001
            velocity(1,i)=0;
        end
    end
    
    velocity2 = zeros(2*npts,1);
    velocity2(1:2:end) = velocity(1,:)' ;
    velocity2(2:2:end) = velocity(2,:)' ;
    velocity=velocity2;
elseif strcmp(problem,'Cylinder')
    for ii=1:npts
        
        if xpts(ii)<0.000001 || xpts(ii)>2.1999
            % inflow and outflow
            velocity(1,ii)= 6/.41/.41*ypts(ii)*(.41-ypts(ii));
            velocity(2,ii)=0;
        end
    end    
    velocity2 = zeros(2*npts,1);
    velocity2(1:2:end) = velocity(1,:)' ;
    velocity2(2:2:end) = velocity(2,:)' ;
    velocity=velocity2;
    
else    
    error('initvel: specified problem invalid')
end
