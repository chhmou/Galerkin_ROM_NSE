
figure


[X,Y]=meshgrid(0:.02:1,0:.02:1);
    
Uvel = griddata(nodeco(:,1),nodeco(:,2),GlobalV(1:2:end,1),X,Y);
Vvel = griddata(nodeco(:,1),nodeco(:,2),GlobalV(2:2:end,1),X,Y);
quiver(X,Y,Uvel,Vvel,3,'k');

axis tight
