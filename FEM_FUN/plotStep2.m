
% datafile.m must contain nodeco, GlobalV, numXcos, numYcos
% nodeco is nx2, where each row is the coords for a node
% GlobalV holds the velocites (both x and y) as a vector. 
% GlobalV = [u1(node1), u2(node1),u1(node2),u2(node2),...]
% numXcos is the number of x coordinates
% numYcos is the number of y coordinates
figure

 sx = [ 0,   0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 4.5, 5.5, 6.0, 6.25, 6.5, 8, 9, 10, 11, 14, 14] ;
 sy = [ 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 0.5, 1.1, 1.1, 0.75, 0.5, 0.5, 0.4, 0.6, 0.8, 0.8, 0.5] ;
% l=[0 0.001 0.01 0.02 0.05 0.1 0.3 0.5 0.7 0.9 1 1.1 1.2];

%load(datafile)
% put the data into the (stupid) plotting format
% start with data as each node getting a row, and each row as [xco,yco,u1,u2]
data = [nodeco(:,1),nodeco(:,2),GlobalV(1:2:end,1),GlobalV(2:2:end,1)];
data1 = sortrows(data,2);
data2 = sortrows(data1,1);


[X,Y]=meshgrid(0:.05:40,0:.05:10);
speed = griddata(nodeco(:,1),nodeco(:,2),sqrt(GlobalV(1:2:end,1).^2 + GlobalV(2:2:end,1).^2 ),X,Y);
contourf(X,Y,speed,linspace(0,1,10))
% shading interp;
 colormap jet;
hold on

 Uvel = griddata(nodeco(:,1),nodeco(:,2),GlobalV(1:2:end,1),X,Y);
 Vvel = griddata(nodeco(:,1),nodeco(:,2),GlobalV(2:2:end,1),X,Y);
 
 h=streamline(X,Y,Uvel,Vvel ,sx,sy) ;
 set(h,'Color','white','LineWidth',2);
 hold on

for yy=0:.01:1
    plot( [5:.01:6],yy*ones(101,1),'k')
end
 
colorbar
set(gca,'FontSize',18)
