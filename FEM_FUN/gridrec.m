function [nodeco, elnode, eldegr, elnbgh, bdynde] =  ...
                                          gridrec(xvec, nxvec, yvec, nyvec)
%
% This function generates a criss-cross triangulation of the rectangular
% region described by the vectors xvec and yvec. The entries with xvec and
% yvec define x and y partitions of the region. The arrays nxvec and nyvec
% describe the number of divisions to generate within each of the partitions.
%

[dimxvec, temp] = size(xvec) ;
[dimyvec, temp] = size(yvec) ;

x1vec = [ ] ; y1vec = [ ] ; nodeco = [ ] ;
% Generate x-coordinates.
for i = 1 : dimxvec-1
    tempv = linspace(xvec(i) , xvec(i+1) , nxvec(i)+1)' ;
    x1vec = [x1vec ; tempv(1:nxvec(i))] ;
end
x1vec = [x1vec ; xvec(dimxvec)] ;

[dimx1vec, temp] = size(x1vec) ;
x2vec(1:dimx1vec-1,1) = 0.5*( x1vec(1:dimx1vec-1) + x1vec(2:dimx1vec) ) ;
[dimx2vec, temp] = size(x2vec) ;

%% Generate y-coordinates.
for i = 1 : dimyvec-1
    tempv = linspace(yvec(i) , yvec(i+1) , nyvec(i)+1)' ;
    y1vec = [y1vec ; tempv(1:nyvec(i))] ;
end
y1vec = [y1vec ; yvec(dimyvec)] ; 

[dimy1vec, temp] = size(y1vec) ;
y2vec(1:dimy1vec-1,1) = 0.5*( y1vec(1:dimy1vec-1) + y1vec(2:dimy1vec) ) ;
[dimy2vec, temp] = size(y2vec) ;

%% Set up up the coordinates of the nodes.
for i = 1 : dimy1vec-1
    nodeco = [nodeco ; x1vec, y1vec(i)* ones(dimx1vec,1); ...
                                x2vec, y2vec(i)* ones(dimx2vec,1)]  ;
end
nodeco = [nodeco ; x1vec, y1vec(dimy1vec)* ones(dimx1vec,1)] ;
[nnds, temp] = size(nodeco) ;

%% Now for the definition of the triangles and their neighbours.
ntri = 0 ;
for j = 1 : dimy1vec-1
    for i = 1 : dimx1vec-1
       ll = (j-1)*(2*dimx1vec - 1) + i ;
       lr = ll + 1 ;
       ul = j*(2*dimx1vec - 1) + i ;
       ur = ul + 1 ;
       cc = ll + dimx1vec ;

       ntri = ntri + 1;
       elnode(ntri, 1) = cc ;
       elnode(ntri, 2) = ll ;
       elnode(ntri, 3) = lr ;
%
       elnbgh(ntri,1) = ntri + 3 ;
       if (j == 1) 
           elnbgh(ntri,2) = 0 ;
       else
           elnbgh(ntri,2) = ntri - 4*(dimx1vec-1) + 2 ;
       end
       elnbgh(ntri,3) = ntri + 1 ;
       
       ntri = ntri + 1;
       elnode(ntri, 1) = cc ;
       elnode(ntri, 2) = lr ;
       elnode(ntri, 3) = ur ;
%
       elnbgh(ntri,1) = ntri - 1 ;
       if (i == dimx1vec-1) 
           elnbgh(ntri,2) = 0 ;
       else
           elnbgh(ntri,2) = ntri + 6 ;
       end
       elnbgh(ntri,3) = ntri + 1 ;

       ntri = ntri + 1;
       elnode(ntri, 1) = cc ;
       elnode(ntri, 2) = ur ;
       elnode(ntri, 3) = ul ;
%
       elnbgh(ntri,1) = ntri - 1 ;
       if (j == dimy1vec-1) 
           elnbgh(ntri,2) = 0 ;
       else
           elnbgh(ntri,2) = ntri + 4*(dimx1vec-1) - 2 ;
       end
       elnbgh(ntri,3) = ntri + 1 ;

       ntri = ntri + 1;
       elnode(ntri, 1) = cc ;
       elnode(ntri, 2) = ul ;
       elnode(ntri, 3) = ll ;
%
       elnbgh(ntri,1) = ntri - 1 ;
       if (i == 1) 
           elnbgh(ntri,2) = 0 ;
       else
           elnbgh(ntri,2) = ntri - 6 ;
       end
       elnbgh(ntri,3) = ntri - 3 ;
       
    end
end

%  Initially we assume that out approximation is piecewise linear
%  on all the triangles.
eldegr = ones(ntri,1) ;

% Now we define the boundary nodes and the associated boundary conditions.
% If the node has Dirichlet conditions on the intervals either side: 0
% Dirichlet condition to the left, Neumann to the right: 1
% Neumann to the left and right: 2
% Neumann to the left and Dirichlet to the right: 3
%
% In the case presented here we have that the region represents a 
% rectangular region given by:
%       [xvec(1), xvec(dimxvec)]x[yvec(1), yvec(dimyvec)],
% with Dirichlet boundary conditions all around.
%
%*** NOTE: The boundary is defined counterclockwise.
nbdy = 0 ;
for i = 1:dimx1vec
    nbdy = nbdy + 1 ;
    bdynde(nbdy,1) = i ;
    bdynde(nbdy,2) = 0 ;
end
for j = 2:dimy1vec
    nbdy = nbdy + 1 ;
    bdynde(nbdy,1) = (j-1)*(2*dimx1vec - 1) + dimx1vec ;
    bdynde(nbdy,2) = 0 ;
end
for i = 1:dimx1vec-1
    nbdy = nbdy + 1 ;
    bdynde(nbdy,1) = nnds - i ;
    bdynde(nbdy,2) = 0 ;
end
for j = 2:dimy1vec-1
    nbdy = nbdy + 1 ;
    bdynde(nbdy,1) = 1 + (dimy1vec - j)*(2*dimx1vec - 1) ;
    bdynde(nbdy,2) = 0 ;
end
