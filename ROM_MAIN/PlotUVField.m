function PlotUVField(veluv,nodeco)

%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
figure
%POD_V= zeros(length(GlobalV),1);
%POD_V(DofV) = POD_all(:,plot_r);
%RomUV = PodAll2d(:,plot_r);
%POD_V = PhiR(:,plot_r);
%POD_V(bdrydof) = bdrydirvals;
l=[0 0.01 0.05 0.1 0.15 .2 .25 .3 .4 .5 0.6 0.7 0.8 1 1.2 1.4 1.7 2];
% 

% put the data into the (stupid) plotting format
% start with data as each node getting a row, and each row as [xco,yco,u1,u2]




%subplot(2,1,1)

[X,Y]=meshgrid(0:.02:2.2,0:.01:.41);
    
% Uvel = griddata(nodeco(:,1),nodeco(:,2),POD_V(1:2:end,1),X,Y);
% Vvel = griddata(nodeco(:,1),nodeco(:,2),POD_V(2:2:end,1),X,Y);
% quiver(X,Y,Uvel,Vvel,1,'k');
% 
% axis tight
% hold on
% 
% % plot the circle
% for theta=0:.01:2*pi
%     xend=.2 + .05*cos(theta);
%     yend=.2 + .05*sin(theta);
%     xplot = linspace(.2,xend,50);
%     yplot = linspace(.2,yend,50);
%     plot( xplot,yplot,'k')
% end
% set(gca,'DataAspectRatio',[1 1 1])
% 
% 
% subplot(2,1,2)

speed = griddata(nodeco(:,1),nodeco(:,2),sqrt(veluv(1:2:end,1).^2 + veluv(2:2:end,1).^2 ),X,Y);
contourf(X,Y,speed,100,'edgecolor','none')
shading flat;
colormap jet;
hold on
 
axis tight

% plot the circle
for theta=0:.01:2*pi
    xend=.2 + .05*cos(theta);
    yend=.2 + .05*sin(theta);
    xplot = linspace(.2,xend,50);
    yplot = linspace(.2,yend,50);
    plot( xplot,yplot,'w')
end
set(gca,'DataAspectRatio',[1 1 1])


title(['Velocity field for the flow past the cylinder '])




