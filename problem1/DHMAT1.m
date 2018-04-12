function [A,B] = DHMAT1( W,Z,DT )

Wd=W*sqrt(1-(Z^2));
a11=exp(-Z*W*DT)*((Z/sqrt(1-Z^2))*sin(Wd*DT)+cos(Wd*DT));
a22=exp(-Z*W*DT)*(cos(Wd*DT)-((Z/sqrt(1-Z^2))*sin(Wd*DT)));
a12=((exp(-Z*W*DT))/Wd)*sin(Wd*DT);
a21=(-W/sqrt(1-Z^2))*(exp(-Z*W*DT)*sin(Wd*DT));
b11=exp(-Z*W*DT)*(((((2*(Z^2)-1)/((W^2)*DT))+(Z/W))*(sin(Wd*DT)/Wd))+(((2*Z/((W^3)*DT))+(1/W^2))*cos(Wd*DT)))-((2*Z)/((W^3)*DT));
b12=-exp(-Z*W*DT)*(((((2*(Z^2)-1)/((W^2)*DT)))*(sin(Wd*DT)/Wd))+(((2*Z/((W^3)*DT)))*cos(Wd*DT)))+((2*Z)/((W^3)*DT))-(1/W^2);
b21=exp(-Z*W*DT)*((((2*(Z^2)-1)/((W^2)*DT))+(Z/W))*(cos(Wd*DT)-(Z/sqrt(1-(Z^2)))*sin(Wd*DT))-(((2*Z)/((W^3)*DT))+(1/W^2))*((Wd*sin(Wd*DT))+(Z*W*cos(Wd*DT))))+(1/((W^2)*DT));
b22=-exp(-Z*W*DT)*((((2*(Z^2)-1)/((W^2)*DT)))*(cos(Wd*DT)-(Z/sqrt(1-(Z^2)))*sin(Wd*DT))-(((2*Z)/((W^3)*DT)))*((Wd*sin(Wd*DT))+(Z*W*cos(Wd*DT))))-(1/((W^2)*DT));
A=[a11 a12; a21 a22];
B=[b11 b12; b21 b22];
 
end
