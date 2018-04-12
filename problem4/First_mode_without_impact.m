L1=0.4572; %Length of column 1, 18 in converted to m
L2=0.4572; %Length of column 2, 18 in converted to m
L3=0.4572; %Length of column 3, 18 in converted to m
L4=0.4064; %Length of column 4, 16 in converted to m

E_col=200E9; %elastic modulus of columns (carbon steel)
R_col=4.7625e-3; %radius of columns, 3/16in converted to m
I_col=pi*(R_col^4)*0.25; %moment of inertia of columns

k1_calc= 4*((12*E_col*I_col)/(L1^3));
k2_calc= 4*((12*E_col*I_col)/(L2^3));
k3_calc= 4*((12*E_col*I_col)/(L3^3));
k4_calc= 4*((12*E_col*I_col)/(L4^3));

filename = 'Free Vibration.xlsx';
tol=0.6;
[gamma_1, gamma_2, gamma_3, gamma_4]=dampening_coeff(filename,tol);

K=zeros(4,4);
K(1,1)=k1+k2; K(1,2)=-k2;
K(2,1)=-k2; K(2,2)=k2+k3; K(2,3)=-k3;
K(3,2)=-k3; K(3,3)=k3+k4; K(3,4)=-k4; 
K(4,3)=-k4; K(4,4)=k4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
density_steel=0.284; %lb/in^3
density_al=0.0975; %lb/in^3

lb_to_kg=0.45359237;

volume_base=0.25*24*24;
volume_base_bar1=(2*1*24)*2;
volume_base_bar2=(2*1*20)*2;

volume_floor=(1/8)*23.5*23.5;
volue_frame=(0.359375*23)*4;

volume_channels=(26.5*0.3*0.8) * 6;
volume_t_section=((2.25+2)*0.1)*22.5;

volume_vertical_bar=12*1*0.1;

diameter_ball=0.0508; %2 in converted to m

m_ball=0.535239; %1.18 lb converted to kg
m_clamp=0.358338; %0.79 lb converted to kg

m_regular_clamp=0.96688; %2.1316 lb converted to kg


m1=(volume_floor+volue_frame)*density_al*lb_to_kg;
m2=(volume_floor+volue_frame)*density_al*lb_to_kg;
m3=((volume_floor+volume_channels+volume_t_section)*density_al*lb_to_kg)+(2*m_clamp);
m4=((volume_floor+volume_vertical_bar)*density_al*lb_to_kg)+m_regular_clamp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M=zeros(4,4);
M(1,1)=m1; M(2,2)=m2; M(3,3)=m3; M(4,4)=m4;


B=[m1; m2; m3; m4];

c1=gamma_1*sqrt(4*k1*m1);
c2=gamma_2*sqrt(4*k2*m2);
c3=gamma_3*sqrt(4*k3*m3);
c4=gamma_4*sqrt(4*k4*m4);

C=zeros(4,4);
C(1,1)=c1+c2; C(1,2)=-c2;
C(2,1)=-c2; C(2,2)=c2+c3; C(2,3)=-c3;
C(3,2)=-c3; C(3,3)=c3+c4; C(3,4)=-c4; 
C(4,3)=-c4; C(4,4)=c4;

filename = 'First Mode Without Impact.xlsx';
tol=1e-3;
[time, base_disp, base_accel, floor_1, floor_2, floor_3, floor_4, amp_array,freq_array] = Extract_Data(filename,tol);
dt=time(2)-time(1);

load_matrix=zeros(length(time),length(freq_array));
super_load=zeros(length(time),1);
for i=1:length(freq_array)
    load_matrix(:,i)=-1.75*amp_array(i)*sin(2*pi*freq_array(i)*time(:))'; %1.5 factor to account for amplitude loss during FFT
    super_load(:)=super_load(:)+ (-1.75*amp_array(i)*sin(2*pi*freq_array(i)*time(:)));
end

x0=zeros(4,1); 
v0=zeros(4,1);
% [x,xdot] = rk4_MDOF_multiple_loads( dt, x0, v0, M ,C, K, B, load_matrix);
[x,xdot] = rk4_MDOF( dt, x0, v0, M ,C, K, B, super_load);

accel_1=zeros(1,length(time));
accel_2=zeros(1,length(time));
accel_3=zeros(1,length(time));
accel_4=zeros(1,length(time));
for i=2:size(xdot,1)
   accel_1(i)= (xdot(i,1)-xdot(i-1,1))/dt;
   accel_2(i)= (xdot(i,2)-xdot(i-1,2))/dt;
   accel_3(i)= (xdot(i,3)-xdot(i-1,3))/dt;
   accel_4(i)= (xdot(i,4)-xdot(i-1,4))/dt;
end

figure
hold on
plot(time,accel_1,'b');
plot(time,floor_1,'r');
grid on
title('floor 1');
xlabel('time');
ylabel('acceleration');
legend('Numerical','Experimental');

figure
hold on
plot(time,accel_2,'b');
plot(time,floor_2,'r');
grid on
title('floor 2');
xlabel('time');
ylabel('acceleration');
legend('Numerical','Experimental');

figure
hold on
plot(time,accel_3,'b');
plot(time,floor_3,'r');
grid on
title('floor 3');
xlabel('time');
ylabel('acceleration');
legend('Numerical','Experimental');

figure
hold on
plot(time,accel_4,'b');
plot(time,floor_4,'r');
grid on
title('floor 4');
xlabel('time');
ylabel('acceleration');
legend('Numerical','Experimental');
