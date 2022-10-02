clc
close all
clear all

syms x1 x2

%Функция Химмельблау
%f1=(x1^2+x2-11)^2+(x1+x2^2-7)^2;
%Функция Розенброка
f1=(1-x1)^2+100*(x2-x1^2)^2;
%Квадратичная функция
%f1=-12*x2+4*x1^2+4*x2^2-4*x1*x2;
fx=inline(f1);
fobj=@(x) fx(x(:,1),x(:,2));

grad=gradient(f1);
G=inline(grad);
gradx=@(x) G(x(:,1),x(:,2));

H1=hessian(f1);
HH=inline(H1);
Hx=@(x) HH(x(:,1),x(:,2));

x0=[-1,-1];
maxiter=2;
tol=1e-4;
iter=0;
S=0;

[x_fun,fval,exitflag,output] = fminsearch(fobj,x0);

%построение графика
plot(x0(1),x0(2),'*k');
hold on
 
%Функция Розенброка
[X1,X2]=meshgrid(-4:0.05:4);Z=100*(X2-X1.^2).^2+(X1-1).^2;[c,h]=contour...
(X1,X2,Z,[5 12 70 120 200 300 450 600 800],'blue');
clabel(c,h,[70 300 600],'FontSize',7);

%Функция Химмельблау
% [X1,X2]=meshgrid(-5:0.03:5);Z=(X1.^2+X2-11).^2+(X1+X2.^2-7).^2;[c,h]=contour...
% (X1,X2,Z,[-0.8 -0.2 0.07 1 3 7 15 28 45 65],'red');
% clabel(c,h,[70 300 600],'FontSize',7);

% xi1 = -5:0.5:5;
% xi2 = -5:0.1:5;
% [X1, X2] = meshgrid(xi1,xi2);
% Ff = -12.*X2+4*X1.^2+4*X2.^2-4.*X1.*X2;
% meshc(X1, X2, Ff); 
% contour(X1, X2, Ff);
% hold on


Gpr=-gradx(x0);
x_print=[x0(1)];
y_print=[x0(2)];


while norm(gradx(x0))>tol %&& iter<maxiter

    Gi=-gradx(x0);
    H=Hx(x0);
    bet=norm(Gi).^2./norm(Gpr).^2;
    S=Gi+bet.*S;
    lam=Gi'*S./(S'*H*S);
    Xnew=x0+lam.*S';
    x0=Xnew;
    Gpr=Gi;
    iter=iter+1;
    x_print=[x_print x0(1)];
    y_print=[y_print x0(2)];
    
end

plot(x_print,y_print,'k');
fprintf('Optimal Solution x=[%f,%f]\n',x0(1),x0(2));
fprintf('Optimal value f(x)=%f\n',fobj(x0));
fprintf('Number of iterations=%f\n',iter);

%Результат встроенной функции
fprintf('Optimal Solution from fminsearch x=[%f,%f]\n',x_fun(1),x_fun(2));
fprintf('Optimal value from fminsearch f(x)=%f\n',fval);
output

