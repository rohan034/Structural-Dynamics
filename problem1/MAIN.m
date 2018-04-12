
 
W=2*pi;
Z=0;
DT=0.1;
x=0;
xdot=0;
Tmax=3;

 
[A,B]=DHMAT1(W,Z,DT);
A;
B;
 
m=((Tmax/DT));
n=m+1;
i=0;
t=-1*DT;
Xold=[x;xdot];
X=zeros(n,2);
I=zeros(n,1)
time=zeros(n,1);
results=zeros(n,4)
 
for k=0:m;
    i=i+1;
    t=t+DT;
    [ F ] = EXCIT(t,DT);
    Xnew=(A*Xold)-(B*F)
    
    time(i,1)=t;
    I(i,1)=i;
    X(i,1)=Xold(1);
    X(i,2)=Xold(2);
    
    Xold(1)=Xnew(1);
    Xold(2)=Xnew(2);
end

X;
time;
results(:,1)=I;
results(:,2)=time;
results(:,3)=X(:,1);
results(:,4)=X(:,2);
 
Results=results';
    
fileid=fopen('Results.txt','w');
fprintf(fileid, 'Computer Project 1 Results \n \n \n');
fprintf(fileid, 'I= %f \n T= %f \n Xi= %f \n Xdoti= %f \n \n ',Results);
fclose(fileid); 

