function [nodeco, elnode, bdynde, bdyedge, nVert, nedge] =  midEDGEgen(nodeco, elnode, bdynde) ;


%%%%%%%%%%%%%%%%%%%%
%% This function introduces a midedge node in the triangulation and
%%  the appropriate boundary conditions for that node.
%% This midedge node is stored in the definition of elnode.
%%


ntri = size(elnode,1) ;
nVert = size(nodeco,1) ;
nnds = nVert ;

edgeset = sparse(nVert,nVert) ;
nedge = 0 ;

vec1 = [2, 3, 1, 2] ;
for itrg = 1 : ntri
   for iedge = 1 : 3
      n1 = elnode(itrg, vec1(iedge)) ;
      n2 = elnode(itrg, vec1(iedge+1)) ;
      nmin = min(n1, n2) ;
      nmax = max(n1, n2) ;
      if (edgeset(nmin, nmax) == 0)
         nedge = nedge + 1 ;
         nnds = nnds + 1 ;
         nodeco(nnds,:) = (nodeco(n1,:) + nodeco(n2,:)) / 2.0 ;
         edgeset(nmin, nmax) = nedge ;
         edgenode = nedge ;
      else
         edgenode = edgeset(nmin, nmax) ;
      end
      elnode(itrg,3+iedge) = edgenode ;
      
   end
end



%%% Now to introduce the bdyedge node list  
%%% Note: (i) We assume a triangle vertex lies at the points along the
%%            the boundary where the type of boundary condition changes.
%%        (ii) For an edge node to be Dirichlet both vertex nodes on
%%            either side must be Dirichlet.
%%        (iii) Since pure Neumann corresponds to imposing no "Boundary
%%            condition" we can get away with the min statement to define
%%            the bounday conditions for the velocity and pressure.

% nbdy = size(bdynde,1) ;
bdyedge = [ ] ;
% for iby = 1:nbdy-1
%    n1 = bdynde(iby, 1) ;
%    n2 = bdynde(iby+1, 1) ;
%    nmin = min(n1, n2) ;
%    nmax = max(n1, n2) ;
%    vbc = min( bdynde(iby, 2) , bdynde(iby+1, 2) ) ;
%    sbc = min( bdynde(iby, 4) , bdynde(iby+1, 4) ) ;
%    bdyedge = [bdyedge ; [edgeset(nmin,nmax), vbc, -1, sbc] ] ;
% end
% n1 = bdynde(1, 1) ;
% n2 = bdynde(nbdy, 1) ;
% nmin = min(n1, n2) ;
% nmax = max(n1, n2) ;
% vbc = min( bdynde(1, 2) , bdynde(nbdy, 2) ) ;
% sbc = min( bdynde(1, 4) , bdynde(nbdy, 4) ) ;
% bdyedge = [bdyedge ; [edgeset(nmin,nmax), vbc, -1, sbc] ] ;

         
         
   
   
