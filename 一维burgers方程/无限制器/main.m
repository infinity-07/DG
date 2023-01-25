clc,clear,close all
num = 40;
time = 0.5/pi;

[xx,U_total,delta_t]=dg_solver(num,time);



U = U_total{end};
yy1 = zeros(1,length(xx));
yy2 = zeros(1,length(xx));
syms xii

for i = 1:length(xx)
    x = xx(i);
    yy1(i) = Compute_U(U,i,x);
    x0 = vpasolve(xii+time*sin(xii)-x);
    yy2(i) = double(sin(x0));
end

plot(xx,yy1,'*')
hold on
plot(xx,yy2,'LineWidth',2)

err = norm(yy1-yy2,inf);
disp(err)





