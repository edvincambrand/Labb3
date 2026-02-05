%% Framåt Euler

function w = framat_euler(f, tspan, y0, n)
    h = (tspan(2)-tspan(1))/n;
    p = size(y0);
    w = zeros(p(1), n+1);
    w(:,1) = y0;

    for i=1:n
        w(:,i+1) = w(:,i) + h*f(w(:,i));
    end
end


%Kallar framåt Euler
tspan = [0,100];
a = 0.5;
y0 = [1-a; 0; 0; sqrt((1+a)/(1-a))];
n = 100000; 
f = @(y) [y(3); y(4); -y(1)/((y(1)^2+y(2)^2)^1.5); -y(2)/((y(1)^2+y(2)^2)^1.5)];

y = framat_euler(f,tspan, y0, n);

%Plottar banorna (q1(t),q2(t))
q1 = y(1,:);
q2 = y(2,:);

plot(q1,q2)

grid on


%% Bakåt Euler 


%Definierar Bakåt Euler
function [t,w] = bakat_euler(F,DF,tspan,y0,n)
    h = (tspan(2)-tspan(1))/n;
    t = tspan(1) + h*(0:n);
    p = size(y0);
    w = zeros(p(1), n+1);
    w(:,1) = y0;
    %Använder Newtons Metod approximera w_n+1 
    for i=1:n
        w(:,i+1) = newton(w(:,i),F,DF,h);
    end
end

%Definierar funktionen F(x)=x-y_i-h*f(y_i), där x = y_i+1     
F = @(x,xb,h) [x(1)-xb(1)-h*x(3);
               x(2)-xb(2)-h*x(4);
               x(3)-xb(3)+h*x(1)/(x(1)^2+x(2)^2)^1.5;
               x(4)-xb(4)+h*x(2)/(x(1)^2+x(2)^2)^1.5];

%Definierar Jacobianen av F
DF = @(x,h) [1, 0, -h, 0;
             0, 1, 0, -h;
            h* ((x(1)^2+x(2)^2) - 3*x(1)^2) / ( x(1)^2+x(2)^2 )^2.5, h* (-3*x(1)*x(2)*( x(1)^2+x(2)^2 )^0.5) / ( x(1)^2+x(2)^2 )^3, 1, 0;
            h* (-3*x(1)*x(2)*( x(1)^2+x(2)^2 )^0.5) / ( x(1)^2+x(2)^2 )^3, h* ((x(1)^2+x(2)^2) - 3*x(2)^2) / ( x(1)^2+x(2)^2 )^2.5, 0, 1];

%Definierar Newtons Metod
function y = newton(wb,F,DF,h)
    s = [1; 1; 1; 1];
    tol = 1e-8;
    x = wb;
    while norm(s) >= tol
        s = - DF(x,h) \ F(x,wb,h);
        x = x + s;
    end

    y = x;
end

%Kallar bakåt Euler
tspan = [0,50];
a = 0.5;
y0 = [1-a; 0; 0; sqrt((1+a)/(1-a))];
n = 200000; 
[t,y] = bakat_euler(F,DF,tspan, y0, n);

%Plottar banorna (q1(t),q2(t))
q1 = y(1,:);
q2 = y(2,:);

plot(q1,q2)

grid on

