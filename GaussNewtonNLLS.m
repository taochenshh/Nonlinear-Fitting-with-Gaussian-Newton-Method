function xc=GaussNewtonNLLS(r,x,x0,epsilon,n)
%Gauss-Newton NonLinear Least Squares
%Input: r is the nonlinear fitting function model
%       x is the coefficients in the fitting function that needs to be
%       determined
%       x0 is the initial values
%       epsilon is the tolerance
%       n is the maximum iteration number
%Output: xc is the determined coefficients
x0=x0(:);
x=x(:);
J=jacobian(r,x);
a=1;
A0=double(subs(J,x,x0));
r0=double(subs(r,x,x0));
v=-(A0'*A0)\(A0'*r0);
x0=x0+v;
while(norm(A0'*r0)>epsilon)
    A0=double(subs(J,x,x0));
    r0=double(subs(r,x,x0));  
    v=-(A0'*A0)\(A0'*r0);
    x0=x0+v;
    a=a+1;
    if a>n break;
    end
    a=a+1;
end
xc=x0;
end

