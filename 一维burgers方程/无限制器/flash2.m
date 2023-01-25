figure
hold on
ylim([-1.1,1.1])

h1 = plot(xx,yy1(1,:),'r*');
h2 = plot(xx,yy2(2,:),'b','LineWidth',2);
pause


for i = 1:n
    delete(h1)
    delete(h2)
    h1 = plot(xx,yy1(i,:),'r*');
    h2 = plot(xx,yy2(i,:),'b','LineWidth',2);
    pause(0.01)
end