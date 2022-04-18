function bdryvals = getBdryVals_sc(boundary_dof,time) 

global problem GlobalV nodeco  nu

boundary_dof2 = boundary_dof;
boundary_dof(1:2:end) = boundary_dof2(1:end/2);
boundary_dof(2:2:end) = boundary_dof2(end/2+1:end);
boundary_dof = boundary_dof';


% First we need the xy coordinates for the boundary_dof
% We assume the boundary_dof are in order, so they come in pairs in 2d
if mod(size(boundary_dof,1),2) == 1
    error('size of boundarydof is not right!')
end

% convert from global dof number to node number
% Assumes if one dof from a node is dirichlet, then so is the other.
xypts= nodeco( boundary_dof(2:2:end) / 2 , :);
xpts = xypts(:,1).';
ypts = xypts(:,2).';
npts = size(xpts,2) ;

% Each point corresponds to two bdryvals (the x and y piece)

bdryvals = 0*boundary_dof;

 if strcmp(problem,'Step')
     for i=1:npts
         dofnum=2*i-1;
         if  abs( ypts(i)*(ypts(i)-10)) < 1e-12
             % top and bottom have u1=u2=0
             bdryvals(dofnum)=0;
             bdryvals(dofnum+1)=0;
         elseif abs( xpts(i)*(40-xpts(i)) ) < 1e-12
             % sides
             bdryvals(dofnum)= ypts(i)*(10-ypts(i))/25;
             bdryvals(dofnum+1)=0;
         elseif (xpts(i)<=6 & xpts(i)>=5 & ypts(i)<=1 )
             % step
             bdryvals(dofnum)=0;
             bdryvals(dofnum+1)=0;
         else
             error(['Boundary point not on boundary: x=' num2str(xpts(i)) ', y=' num2str(ypts(i)) ])
         end
     end
 elseif strcmp(problem,'Cylinder')
 
     % Cylinder problem
     % if on inlet or outlet, then need nonzero bc    
     for ii=1:npts
         if abs (xpts(ii) * (xpts(ii)-2.2))<1e-10
             bdryvals(2*ii-1) = (6/(.41^2)) * sin(pi*time/8)*ypts(ii)*(.41-ypts(ii));
         end
     end
 
 elseif strcmp(problem,'trig')
     bdryvals(1:2:end) = (1 + 0.01 * time) * sin(ypts);
     bdryvals(2:2:end) = (1 + 0.01 * time) * cos(xpts);
     
 elseif strcmp(problem,'Stokestrig')
     bdryvals(1:2:end) = cos(2*pi*ypts);
     bdryvals(2:2:end) = sin(2*pi*xpts);
     
 elseif strcmp(problem,'cavity')
     for i=1:npts
         dofnum=2*i-1;
         if  ypts(i)>.9999
             bdryvals(dofnum)=1;
             bdryvals(dofnum+1)=0;
 
         end
     end
 else
     
     error('getbdryvals: invalid problem specified')
%     %% Chorin problem
%     n=2;
%     bdryvals(1:2:end) = -cos(n*pi*xpts).*sin(n*pi*ypts) * exp(-2* n^2 * pi^2 * t * nu);
%     bdryvals(2:2:end) =  sin(n*pi*xpts).*cos(n*pi*ypts) * exp(-2* n^2 * pi^2 * t * nu);
 end
% 

bdryvals2 = bdryvals;
bdryvals(1:end/2) = bdryvals2(1:2:end);
bdryvals(end/2+1:end) = bdryvals2(2:2:end);
