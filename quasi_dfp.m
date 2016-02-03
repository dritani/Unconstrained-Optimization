THEX = 0;
THEY = 0;
CONV = 0;
% These are only here for quick reference
f = (x - 1)^2 + 100*(- x^2 + y)^2;
grad = [2*x - 400*x*(- x^2 + y) - 2 ; - 200*x^2 + 200*y];
p_k =  [400*x*(- x^2 + y) - 2*x + 2;  200*x^2 - 200*y];
the_hess = [ 1200*x^2 - 400*y + 2, -400*x;-400*x,    200];

k=0;
limit = 100000;
tol = 0.01;
alpha = 0.001;
theC = 10^(-4);
rho = 0.9;

H_new = eye(2);
x = 3;
y = 4;

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
    H_old = H_new;
    
    temp = -H_old * [2*x - 400*x*(- x^2 + y) - 2 ; - 200*x^2 + 200*y];
    x = x_old + alpha * temp(1);
    y = y_old + alpha * temp(2);
    

    deltaX = alpha * temp;
    %deltaG = grad(x,y) - grad(x_old,y_old);
    deltaG = [2*x - 400*x*(- x^2 + y) - 2 ; - 200*x^2 + 200*y] - ...
        [2*x_old - 400*x_old*(- x_old^2 + y_old) - 2 ; - 200*x_old^2 + 200*y_old];
    
    H_new = H_old + ((deltaX * deltaX')/(deltaX' * deltaG)) - ...
        (H_old * deltaG * deltaG' * H_old )/(deltaG' * H_old * deltaG);


end

x
y
k
epsilon

semilogy(CONV)
xlabel('Number of iterations')
ylabel('Norm of \nablaf(x,y)')
title('Figure 7 - Quasi-Newton DFP - Convergence of gradient')
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
title('Figure 8 - Quasi-Newton DFP - Contour Plot and Algorithm Path')

dlmwrite('qu_CONV.txt',CONV);