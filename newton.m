THEX = 0;
THEY = 0;
CONV = 0;
syms x y;
syms p_k;
syms grad;
f(x,y) = (1-x)^2+100*((y-(x^2))^2);
grad(x,y) = gradient(f);
hess(x,y) = hessian(f);
p_k(x,y) = -gradient(f);

i=0;
limit = 10;
tol = 0.01;
alpha = 0.0005;
theC = 10^(-4);
rho = 0.9;

x_old = 3;
y_old = 4;
x = 3;
y = 4;

while(1)
    epsilon = norm(double(grad(x,y)));
    
    i=i+1;

    CONV(i) = epsilon;
    THEX(i) = x;
    THEY(i) = y;
    
    if epsilon < tol
        break;
    end  
    
    if (i==limit)
        break;
    end 
    
    x_old = x;
    y_old = y;
    
    delta = -inv(hess(x,y))*grad(x,y);   
    delta = double(delta);
    x = x_old + delta(1);
    y = y_old + delta(2);
end

x
y
i
epsilon

semilogy(CONV)
xlabel('Number of iterations')
ylabel('Norm of \nablaf(x,y)')
title('Figure 5 - Newton''s Method - Convergence of gradient')
grid on
size = 10;
a = linspace(-size,size);
b = linspace(-size,size);
[A,B] = meshgrid(a,b);
C = (1-A).^2+100*((B-(A.^2)).^2);
levels = 100:100:100;
figure
hold on
contour(A,B,C,200)
plot(THEX,THEY)
plot(3,4,'g*')
plot(x,y,'r*')
legend('\nablaf(x,y)','Path','Start','End')
xlabel('X')
ylabel('Y')
title('Figure 6 - Newton''s Method - Contour Plot and Algorithm Path')

dlmwrite('ne_CONV.txt',CONV);