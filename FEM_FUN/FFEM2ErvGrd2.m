function [nodeco, elnode, bdynde] = FFEM2ErvGrd2(ffemgridfile) 
%
%  Revised 09/14/09 --- replaces FFEM2ErvGrd
%  This revision corrects an error in mapping the boundary from FreeFEM to
%  the format needed for my codes
%
% This function takes the data from the grid file generated using FreeFEM
% and converts it into a format needed for our routines. 
% Following this file we need to run it through a file to setup the
% (i) boundary conditions for the nodes on the boundary
% (ii) a mid-edge generator 
%

fin = fopen(ffemgridfile,'r') ;

% Determine nVert, ntri, nbdy
buffer = fgetl(fin) ;
n1 = str2num(buffer) ;
nVert = n1(1) ;
ntri = n1(2) ;
nbdy = n1(3) ;
%[nVert, ntri, nbdy] = [n1(1), n1(2), n1(3)] ;

nodeco = zeros(nVert,2) ;
elnode = zeros(ntri,3) ;
bdynde = zeros(nbdy,1) ;

for ii = 1:nVert
    buffer = fgetl(fin) ;
    vert = str2num(buffer) ;
    nodeco(ii,1) = vert(1) ;
    nodeco(ii,2) = vert(2) ;
end

for ii = 1:ntri
    buffer = fgetl(fin) ;
    n1 = str2num(buffer) ;
    elnode(ii,1) = n1(1) ;
    elnode(ii,2) = n1(2) ;
    elnode(ii,3) = n1(3) ;
end


%%%%%
%%%%%
%% This next part is complicated by the fact that FreeFEM lists boundary segements
%% pieces and does not list them in a continuous fashion around the boundary.

bdymat = zeros(nbdy , 3);  
nedgept = zeros(nbdy, 1) ;


buffer = fgetl(fin) ;
n1 = str2num(buffer) ;
bdymat(1,:) = n1 ;
nedgept(1) = 1 ;
edgeind = n1(3) ;

edgetot = 1 ;

for ii = 2:nbdy
    buffer = fgetl(fin) ;
    n1 = str2num(buffer) ;
    bdymat(ii,:) = n1 ;
%
    if (n1(3) ~= edgeind)
%%      new edge detected
        nedgept(2*edgetot) = ii - 1 ;
        edgetot = edgetot + 1 ;
        nedgept(2*edgetot - 1) = ii ;
%
        edgeind = n1(3) ;
    end
end
nedgept(2*edgetot) = nbdy ;

%%
%% First we order the boundary segments
bdyord = [ ] ;
for ii = 1:edgetot
    ifnd = 0 ;
    ij = 1 ;
    while ( (ij <= edgetot) && (ifnd == 0) )
       if ( ii == bdymat(nedgept(2*ij - 1), 3) )
          bdyord = [bdyord ; bdymat(nedgept(2*ij - 1) : nedgept(2*ij) , :)] ;
          ifnd = 1 ;
       end
       ij = ij + 1 ;
    end
end

bdymat = bdyord ;

%% At each boundary segment we must check if the segment is listed in
%% the correct direction or reverse direction

% Handling the first edge is easy
bdynde(1:nedgept(2)) = bdymat(1:nedgept(2), 1) ;

%% The next node along the boundary is ndenxt


% for ii = 2:edgetot
%    ndenxt = bdymat(nedgept(2*(ii -1)), 2) ;
% 
%    if ( bdymat(nedgept(2*ii - 1), 1) == ndenxt )
%       %% boundary segment is listed consistently
%       bdynde(nedgept(2*ii - 1):nedgept(2*ii)) = bdymat(nedgept(2*ii - 1):nedgept(2*ii), 1) ;
%    else
%       %% segment is listed backwards
%       bdynde(nedgept(2*ii - 1):nedgept(2*ii)) = bdymat(nedgept(2*ii):-1:nedgept(2*ii - 1), 1) ;
%    end
% 
% end
    

st = fclose(fin) ;
