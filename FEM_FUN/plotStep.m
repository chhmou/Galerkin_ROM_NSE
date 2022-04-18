function plotStep(datafile)
% datafile.m must contain nodeco, GlobalV, numXcos, numYcos
% nodeco is nx2, where each row is the coords for a node
% GlobalV holds the velocites (both x and y) as a vector. 
% GlobalV = [u1(node1), u2(node1),u1(node2),u2(node2),...]
% numXcos is the number of x coordinates
% numYcos is the number of y coordinates

 sx = [ 0,   0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 4.5, 5.5, 6.0, 6.25, 6.5, 8, 9, 10, 11, 14, 14] ;
 sy = [ 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 0.5, 1.1, 1.1, 0.75, 0.5, 0.5, 0.4, 0.6, 0.8, 0.8, 0.5] ;
 l=[0 0.001 0.01 0.02 0.05 0.1 0.3 0.5 0.7 0.9 1 1.1 1.2];

load(datafile)
% put the data into the (stupid) plotting format
% start with data as each node getting a row, and each row as [xco,yco,u1,u2]
data = [nodeco(:,1),nodeco(:,2),GlobalV(1:2:end,1),GlobalV(2:2:end,1)];
data1 = sortrows(data,2);
data2 = sortrows(data1,1);


% need to know how many unique x and y coords
xcos = [nodeco(1,1)];
ycos = [nodeco(1,2)];
for ii=2:size(nodeco,1);
    x=nodeco(ii,1);
    y=nodeco(ii,2);
    newX=1;
    newY=1;
    for jj=1:size(xcos,1)
        if abs(x - xcos(jj))<1e-13
            newX=0;
        end
    end
    for jj=1:size(ycos,1)
        if abs(y - ycos(jj))<1e-13
            newY=0;
        end
    end
    if newX==1
        xcos=[xcos;x];
    end
    if newY==1
        ycos=[ycos;y];
    end
end
numXcos = size(xcos,1);
numYcos = size(ycos,1);



xco = reshape(data2(:,1),numYcos,numXcos);
yco = reshape(data2(:,2),numYcos,numXcos);
uvel = reshape(data2(:,3),numYcos,numXcos);
vvel = reshape(data2(:,4),numYcos,numXcos);

speed =  (uvel.^2+vvel.^2).^(1/2);
  
contourf(xco,yco,speed,l);
shading flat;
colormap jet;
colorbar;
hold on;
%h=quiver( xcoord, ycoord, uvel,vvel );
h=streamline(xco,yco,uvel,vvel,sx,sy) ;
set(h,'Color','white');

for yy=0:.01:1
    plot( [5:.01:6],yy*ones(101,1),'k')
end
 
