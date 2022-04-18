figure('units','normalized','outerposition',[0 0 1 1])
for i=1:size(elnode,1)
    nodes = elnode(i,1:3);
    nodes(4)=nodes(1);
    for j=1:3
        n1 = nodes(j);
        n2 = nodes(j+1);
        xy1 = nodeco(n1,1:2);
        xy2 = nodeco(n2,1:2);
        xplot=linspace(xy1(1),xy2(1),100);
        yplot=linspace(xy1(2),xy2(2),100);
        plot(xplot,yplot,'k')
        hold on
    end
    
end
xlim([min(nodeco(:,1)) max(nodeco(:,1))])
ylim([min(nodeco(:,2)) max(nodeco(:,2))])
set(gca,'DataAspectRatio',[1 1 1])
title('2D FEM mesh for the flow past the cylinder ','FontSize',28)
