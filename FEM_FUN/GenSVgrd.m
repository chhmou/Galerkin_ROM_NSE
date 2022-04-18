function zz = GenSVgrd(ffemgridfile) 
%
% This function takes the data from the grid file ffemgridfile generated using FreeFEM
% and converts it into a format needed for Scott-Vogelius elements. 
% The new triangulation is output to the file SVgrid.msh
% 
%
zz=1;
fin = fopen(ffemgridfile,'r') ;

% Determine nVert, ntri, nbdy
buffer = fgetl(fin) ;
n1 = str2num(buffer) ;
nVert = n1(1) ;
ntri = n1(2) ;
nbdy = n1(3) ;
%[nVert, ntri, nbdy] = [n1(1), n1(2), n1(3)] ;

nodes = zeros(nVert + ntri,3) ;
elements = zeros(ntri,4) ;
bdysegs = zeros(nbdy,3) ;

for ii = 1:nVert
    buffer = fgetl(fin) ;
    vert = str2num(buffer) ;
    nodes(ii,1) = vert(1) ;
    nodes(ii,2) = vert(2) ;
    nodes(ii,3) = vert(3) ;
end

for ii = 1:ntri
    buffer = fgetl(fin) ;
    n1 = str2num(buffer) ;
    nn = nVert + ii ;
    nodes(nn, 1) = ( nodes(n1(1), 1) + nodes(n1(2), 1) + nodes(n1(3), 1) ) / 3.0 ;
    nodes(nn, 2) = ( nodes(n1(1), 2) + nodes(n1(2), 2) + nodes(n1(3), 2) ) / 3.0 ;
    nodes(nn, 3) = 0 ;
    nt = 3*(ii - 1) + 1 ;
    elements(nt, 1:4) = [n1(1) n1(2) nn n1(4)] ;
    elements(nt + 1, 1:4) = [n1(2) n1(3) nn n1(4)] ;
    elements(nt + 2, 1:4) = [n1(3) n1(1) nn n1(4)] ;
end

for ii = 1:nbdy
    buffer = fgetl(fin) ;
    n1 = str2num(buffer) ;
    bdysegs(ii,1) = n1(1) ;
    bdysegs(ii,2) = n1(2) ;
    bdysegs(ii,3) = n1(3) ;
end

fclose(fin) ;



fid = fopen('SVgrid.msh', 'w');
%% Write out the number of nodes, elements and boundary segments
fprintf(fid, '%d %d %d\n', [nVert+ntri 3*ntri nbdy]);
fprintf(fid, '%d %d %d\n', nodes');
fprintf(fid, '%d %d %d %d\n', elements');
fprintf(fid, '%d %d %d\n', bdysegs');

fclose(fid);
