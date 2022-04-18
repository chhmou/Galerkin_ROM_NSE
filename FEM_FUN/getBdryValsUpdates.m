function bdryvals = getBdryValsUpdates(boundary_dof,time) 

global problem nodeco t dt


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

if strcmp(problem,'YosidaTest2')
    bdryvals(1:2:end) = cos(ypts) * ( (t+exp(t)) - (t-dt+exp(t-dt)) );
    bdryvals(2:2:end) = sin(xpts) * ( (t+exp(t)) - (t-dt+exp(t-dt)) ); 
elseif strcmp(problem,'Step')
    % zero
else
     error('getbdryvals: invalid problem specified')
end
