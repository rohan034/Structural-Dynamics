function [time, base_disp, base_accel, floor_1, floor_2, floor_3, floor_4] = Extract_Data_Dampening(filename, tol)

sheet1 = 1;

raw_data=xlsread(filename,sheet1);

coefficients=[0.3999, 0.3999, 0.3986, 0.3987, 1.017, 1.2407];

time=zeros(1,size(raw_data,1));
base_disp=zeros(1,size(raw_data,1));
base_accel=zeros(1,size(raw_data,1));
floor_1=zeros(1,size(raw_data,1));
floor_2=zeros(1,size(raw_data,1));
floor_3=zeros(1,size(raw_data,1));
floor_4=zeros(1,size(raw_data,1));

g=9.8;
counter=1;
for i=1:size(raw_data,1)
    time(counter)=raw_data(i,1);
    floor_4(counter)=(raw_data(i,2)/coefficients(1))*g;
    floor_3(counter)=(raw_data(i,3)/coefficients(2))*g;
    floor_2(counter)=(raw_data(i,4)/coefficients(3))*g;
    floor_1(counter)=(raw_data(i,5)/coefficients(4))*g;
    base_accel(counter)=(raw_data(i,6)/coefficients(5))*g;
    base_disp(counter)=raw_data(i,7)/coefficients(6);
    counter=counter+1;
end


