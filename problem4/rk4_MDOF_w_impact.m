function [x,xdot] = rk4_MDOF_w_impact( dt, x0, v0, M ,C, K, B, f,Y, Ydot,gap,e,u)

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

%Store ball location
x_ball_store=x0(5:8);

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
    
    %UPDATE BALL POSITION AND VELOCITY HERE
    xdot_vector(5:8)=Ydot-xdot_vector(3);
    x_ball_delta=xdot_vector(5:8)*dt;
    x_vector(5:8)=x_ball_store+x_ball_delta;
%     x_vector(5:8)=x_vector(5:8)-x_vector(3);
    x_ball_store=x_vector(5:8);
    

    %ADD IMPACT CODITION HERE
    for i=5:8
        if abs(x_vector(i))>=gap
            
            %EQN 9 and 10 of "Effect of Imapct Damper on SDOF" paper
            Ydot_old=Ydot;
            Ydot(i-4)=(((1+e)/(1+u))*xdot_vector(3)) + (((u-e)/(1+u))*Ydot(i-4));
            xdot_vector(3)=((e*(1+u)/(1+e))*Ydot_old(i-4)) + (((1-(u*e))/(1+e))*Ydot(i-4));
            
            fprintf('Impact at ball %d. ',i);
            if x_vector(i)<0
                disp('Impact occured at inner stop');
            else
                disp('Impact occured at outer stop');
            end
            
            x_vector(i)=sign(x_vector(i))*gap;
            x_ball_store=x_vector(5:8);
        end
        
    end
    x(k+1,:)=x_vector';
    xdot(k+1,:)=xdot_vector';
    X=[x_vector;xdot_vector];
end
end