%% This is the main program
clear all;
clc;
close all;

n=3; % Number of DOF's
k=3; % Number of modes that are of interest
m1=1; m2=1; m3=2;          % Mass information
k1=1; k2=1; k3=2;          % Stiffness information
c1=0; c2=0; c3=0;          % Damping information
M =[m1 0 0; 0 m2 0;0 0 m3];
K =[k1+k2 -k2 0; -k2 k2+k3 -k3;0 -k3 k3];

x0=[3;1;-1];
[n,M,K,K_inv,k,w,X]=powermethod(n,M,K,k,x0);
x0=[1;1;1];
[n,M,K,K_inv,k,w,X]=powermethod(n,M,K,k,x0);
x0=[-2;1;-1];
[n,M,K,K_inv,k,w,X]=powermethod(n,M,K,k,x0);

disp('Natural frequency for 1st mode (rad/s)');disp(w(1));
disp('Mode shape for 1st mode');disp(X(:,1));


disp('Natural frequencies for 2nd mode (rad/s)');disp(w(2));
disp('Mode shape for 2nd mode');disp(X(:,2));

disp('Natural frequencies for 3rd mode (rad/s)');disp(w(3));
disp('Mode shape for 3rd mode');disp(X(:,3));

