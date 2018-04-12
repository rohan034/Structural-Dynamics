function [n,M,K,K_inv,k,w,X]=powermethod(n,M,K,k,x0)
% Number of DOF's: n
% Mass matrix: M
% Stiffness matrix: K
% Nubmer of modes: k
% Initial guess: x0
if k>2*n
    disp('Error: Number of modes required needs to be smaller');
    return;
end
K_inv=inv(K);
tol=0.0000000001; %This is the tolerance in the iterations.

x=x0./norm(x0);
D=inv(K)*M;
num=1;

for i=1:k
    while abs(abs(x'*(D*x)/(sqrt(x'*x))/(sqrt((D*x)'*(D*x))))-1)>tol 
        norm_value(num,i)=abs(abs(x'*(D*x)/(sqrt(x'*x))/(sqrt((D*x)'*(D*x))))-1);
        x=D*x;
        num=num+1;
        if num>999999
           disp('Warning: Too many iterations');
           return;
        end 
    end
    x_record=x/sqrt(x'*M*x);
    x_record'*M*x_record;
    X(:,i)=x_record;
    lambda=((inv(K)*M*x)'*x)/(x'*x);
    w(i)=sqrt(1/lambda);
    D=D-lambda*x_record*x_record'*M;
    x=x0;
    num=1;
    
end

% Plot the convergence rate.
figure();
plot(norm_value(:,1));hold on;
plot(norm_value(:,2));hold on;
plot(norm_value(:,3));hold off;
xlabel('number of iterations');
ylabel('error (as defined in this program)');
legend('mode 1','mode 2','mode 3');
title(strcat('initial vector x0=',num2str(x0(1)),';',num2str(x0(2)),';',num2str(x0(3))));
end

