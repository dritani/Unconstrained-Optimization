THEX = 0;
THEY = 0;
CONV = 0;
x = 1;
y = 1;
% These are only here for quick reference
f = (x - 1)^2 + 100*(- x^2 + y)^2;
grad = [2*x - 400*x*(- x^2 + y) - 2 ; - 200*x^2 + 200*y];
p_k_formula =  [400*x*(- x^2 + y) - 2*x + 2;  200*x^2 - 200*y];
hess = [ 1200*x^2 - 400*y + 2, -400*x; -400*x,    200];

k=0;
limit = 100000;
tol = 0.01;
alpha = 0.0005;
theC = 10^(-4);
rho = 0.9;

x = 3;
y = 4;
    p_k = [400*x*(- x^2 + y) - 2*x + 2;  200*x^2 - 200*y];

while(1)
    epsilon = norm([2*x - 400*x*(- x^2 + y) - 2 ; - 200*x^2 + 200*y]);
    
    k=k+1;

    CONV(k) = epsilon;
    THEX(k) = x;
    THEY(k) = y;
    
    if (epsilon < tol)
        break;
    end  
    
    if (k==limit)
        break;
    end 
   
    x_old = x;
    y_old = y;

    x = x_old + alpha * p_k(1);
    y = y_old + alpha * p_k(2);
    p_k = [400*x*(- x^2 + y) - 2*x + 2;  200*x^2 - 200*y];

end

x
y
k
epsilon

semilogy(CONV)
xlabel('Number of iterations')
ylabel('Norm of \nablaf(x,y)')
title('Figure 1 - Steepest Descent - Convergence of gradient')
grid on
size = 5;
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
title('Figure 2 - Steepest Descent - Contour Plot and Algorithm Path')

dlmwrite('st_CONV.txt',CONV);