function [gamma_1, gamma_2, gamma_3, gamma_4]=dampening_coeff(filename,tol)
L1=0.4572; %Length of column 1, 18 in converted to m
L2=0.4572; %Length of column 2, 18 in converted to m
L3=0.4572; %Length of column 3, 18 in converted to m
L4=0.4064; %Length of column 4, 16 in converted to m

E_col=200E9; %elastic modulus of columns (carbon steel)
R_col=4.7625e-3; %radius of columns, 3/16in converted to m
I_col=pi*(R_col^4)*0.25; %moment of inertia of columns

k1= 4*((12*E_col*I_col)/(L1^3));
k2= 4*((12*E_col*I_col)/(L2^3));
k3= 4*((12*E_col*I_col)/(L3^3));
k4= 4*((12*E_col*I_col)/(L4^3));

density_steel=0.284; %lb/in^3
density_al=0.0975; %lb/in^3

lb_to_kg=0.45359237;

volume_rods1=(pi*(R_col^2)*L1)*4;
volume_rods2=(pi*(R_col^2)*L2)*4;
volume_rods3=(pi*(R_col^2)*L3)*4;
volume_rods4=(pi*(R_col^2)*L4)*4;

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
m3=((volume_floor+volume_channels+volume_t_section)*density_al*lb_to_kg)+(2*m_clamp);% +(4*m_ball);
m4=((volume_floor+volume_vertical_bar)*density_al*lb_to_kg)+m_regular_clamp;

w1=sqrt(k1/(m1));
w2=sqrt(k2/(m2));
w3=sqrt(k3/(m3));
w4=sqrt(k4/m4);

[time, base_disp, base_accel, floor_1, floor_2, floor_3, floor_4] = Extract_Data_Dampening(filename,tol);

floor_1=floor_1+abs(floor_1(1));
floor_1 = smooth(floor_1);

floor_1_disp=floor_1;
[pks1,locs1] = findpeaks(floor_1_disp,'MinPeakHeight',tol);
index=0;
counter=1;
for i=1:(length(pks1)-1)
    if time(locs1(i+1))-time(locs1(i))<0.3
        if pks1(i)>pks1(i+1)
            index(counter)=i+1;
        else
            index(counter)=i;
        end
        counter=counter+1;
    end
end

if index ~=0
    pks1(index)=[];
    locs1(index)=[];
end


log_dec_1=(log(pks1(10)/pks1(end)))/(length(pks1)-10);
gamma_1=sqrt((log_dec_1^2)/((4*pi*pi)+(log_dec_1^2)));
damp_trend_1=exp(-gamma_1*w1*(time));

% figure
% hold on
% plot(time,floor_1);
% plot(time(locs1),pks1,'ro');
% grid on
% title('Floor 1 Accel');
% figure
% hold on
% plot(time,damp_trend_1,'r--');
% grid on
% title('Logarithmic Decay - Floor 1');

floor_2=floor_2-abs(floor_2(1));
floor_2 = smooth(floor_2);

floor_2_disp=floor_2;
[pks2,locs2] = findpeaks(floor_2_disp,'MinPeakHeight',tol);
index=0;
counter=1;
for i=1:(length(pks2)-1)
    if time(locs2(i+1))-time(locs2(i))<0.2
        if pks2(i)>pks2(i+1)
            index(counter)=i+1;
        else
            index(counter)=i;
        end
        counter=counter+1;
    end
end

if index ~=0
    pks2(index)=[];
    locs2(index)=[];
end

log_dec_2=(log(pks2(6)/pks2(end)))/(length(pks2)-6);
gamma_2=sqrt((log_dec_2^2)/((4*pi*pi)+(log_dec_2^2)));
damp_trend_2=exp(-gamma_2*w2*(time));

% figure
% hold on
% plot(time,floor_2);
% plot(time(locs2),pks2,'ro');
% grid on
% title('Floor 2 Accel');
% figure
% hold on
% plot(time,damp_trend_2,'r--');
% grid on
% title('Logarithmic Decay - Floor 2');


floor_3=floor_3+abs(floor_3(1));
floor_3 = smooth(floor_3);

floor_3_disp=floor_3;
[pks3,locs3] = findpeaks(floor_3_disp,'MinPeakHeight',tol);

log_dec_3=(log(pks3(1)/pks3(end)))/(length(pks3)-1);
gamma_3=sqrt((log_dec_3^2)/((4*pi*pi)+(log_dec_3^2)));
damp_trend_3=exp(-gamma_3*w3*(time));

% figure
% hold on
% plot(time,floor_3);
% plot(time(locs3),pks3,'ro');
% grid on
% title('Floor 13 Accel');
% figure
% hold on
% plot(time,damp_trend_3,'r--');
% grid on
% title('Logarithmic Decay - Floor 3');

floor_4=floor_4+abs(floor_4(1));
floor_4 = smooth(floor_4);

floor_4_disp=floor_4;
[pks4,locs4] = findpeaks(floor_4_disp,'MinPeakHeight',tol);

index=0;
counter=1;
for i=1:(length(pks4)-1)
    if time(locs4(i+1))-time(locs4(i))<0.2
        if pks4(i)>pks4(i+1)
            index(counter)=i+1;
        else
            index(counter)=i;
        end
        counter=counter+1;
    end
end

if index ~=0
    pks4(index)=[];
    locs4(index)=[];
end

log_dec_4=(log(pks4(1)/pks4(end)))/(length(pks4)-1);
gamma_4=sqrt((log_dec_4^2)/((4*pi*pi)+(log_dec_4^2)));
damp_trend_4=exp(-gamma_4*w4*(time));

% figure
% hold on
% plot(time,floor_4);
% plot(time(locs4),pks4,'ro');
% grid on
% title('Floor 4 Accel');
% figure
% hold on
% plot(time,damp_trend_4,'r--');
% grid on
% title('Logarithmic Decay - Floor 4');
