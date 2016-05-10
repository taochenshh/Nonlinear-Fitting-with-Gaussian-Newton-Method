clear all;clc;close all;
syms c1 c2 t y
Y=c1*t*exp(c2*t)-y;
x=1:8;
x=x';
z=[8.0 12.3 15.5 16.8 17.1 15.8 15.2 14.0];
z=z';
r=subs(Y,{t,y},{x,z});
c=[c1;c2];
epsilon=0.5e-8;
c0=[5;0];
n=30;
c_root=GaussNewtonNLLS(r,c,c0,epsilon,n);
display('The coefficients found by using Gauss-Newton method are:');
display(['c=[c1;c2]=' mat2str(c_root)]);
f=c_root(1).*x.*exp(c_root(2).*x);
error=f-z;
RMSE=norm(error,2);
display(['RMSE=' num2str(RMSE)]);
plot(x,z,'o','MarkerFaceColor','g','MarkerSize',8);
hold on;
t=0:0.1:10;
g=c_root(1).*t.*exp(c_root(2).*t);
plot(t,g,'b');
title('Measured level of drug in blood and its fitting curve');
xlabel('hour');
ylabel('concentration (ng/ml)');



