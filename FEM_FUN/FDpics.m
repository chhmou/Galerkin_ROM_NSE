figure
x = [pi/2 - .15,pi/2-.1, pi/2-.05, pi/2, pi/2+0.05 ]
y = sin(x);
plot(x,y,'ko','LineWidth',3)
set(gca,'Ytick',[])
set(gca,'Xtick',[])

xlabel('              x_{i-2}          x_{i-1}         x_{i}           x_{i+1}         x_{i+2}  ','FontSize',16)
set(gca,'FontSize',16)

figure
plot([x(2),x(3)],[y(2),y(3)],'b--','LineWidth',2)
hold on
plot([x(3),x(4)],[y(3),y(4)],'r-.','LineWidth',2)
hold on
plot(x,y,'ko','LineWidth',2)
set(gca,'Ytick',[])
set(gca,'Xtick',[])

xlabel('              x_{i-2}          x_{i-1}         x_{i}           x_{i+1}         x_{i+2}  ','FontSize',16)
legend('slope=backward difference','slope=forward difference')
set(gca,'FontSize',16)
hold on

figure
plot([x(2),x(4)],[y(2),y(4)],'b--','LineWidth',2)
hold on
plot(x,y,'ko','LineWidth',2)
set(gca,'Ytick',[])
set(gca,'Xtick',[])

xlabel('              x_{i-2}          x_{i-1}         x_{i}           x_{i+1}         x_{i+2}  ','FontSize',16)
legend('slope=centered difference')
set(gca,'FontSize',16)
hold on
