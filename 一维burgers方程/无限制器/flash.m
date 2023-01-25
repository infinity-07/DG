n = length(U_total);
u0 = @(x) sin(x);   % 初值
yy1 = zeros(n,length(xx));
yy2 = zeros(n,length(xx));
now = 0;

delta_t = 0.1*2*pi/num;     % 时间步长

h = waitbar(0,'Please wait...');
for i = 1:n
    U = U_total{i};
    if now+delta_t>time
        delta_t = time-now;
    end
    for j = 1:length(xx)
        x = xx(j);
        yy1(i,j) = Compute_U(U,j,x);
        waitbar(i/n,h,'Processing...')
        syms xii
        eqa = xii + now*u0(xii) - x;
        yy2(i,j) = u0(vpasolve(eqa));
    end
    now = now+delta_t;
end
close(h)


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
    pause(0.1)
end