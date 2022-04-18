for i=1:size(bdryelts,1)
    nodes = elnode(bdryelts(i,1),1:3);
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