function [x,xdot] = rk4_MDOF( dt, x0, v0, M ,C, K, B, f)

DOF=size(K,1);

%Constructing A matrix for system of equations
A = zeros(2*DOF,2*DOF);
A(1:DOF,(DOF+1):(2*DOF))=eye(DOF);
A((DOF+1):(2*DOF),1:DOF)=-inv(M)*K;
A((DOF+1):(2*DOF),(DOF+1):(2*DOF))=-inv(M)*C;

N=size(f,1);

%Constructing X vector for system of equations
X = [x0;v0];
x = zeros(N,DOF);
xdot = zeros(N,DOF);
x(1,:)=x0';
xdot(1,:)=v0';

%Solve system of equations using RK-4 method
for k = 1 : N-1
    force1=f(k,:)';
    force2_3=((f(k+1,:)+f(k,:))/2)';
    force4=f(k+1,:)';
    
    f1=B*force1;
    f2_3=B*force2_3;
    f4=B*force4;
    
    g1 = EQNS(X , A, M, DOF, f1);
    g2 = EQNS(X+dt/2*g1, A, M, DOF, f2_3);
    g3 = EQNS(X+dt/2*g2, A, M, DOF, f2_3);
    g4 = EQNS(X+dt*g3, A, M, DOF, f4);
    X = X + dt/6*(g1 + 2*g2 + 2*g3 + g4 );
    
    x_vector=X(1:DOF);
    xdot_vector=X((DOF+1):(2*DOF));
    x(k+1,:)=x_vector';
    xdot(k+1,:)=xdot_vector';
    X=[x_vector;xdot_vector];
end

end

