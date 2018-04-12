function [time, base_disp, base_accel, floor_1, floor_2, floor_3, floor_4, amp_array,freq_array] = Extract_Data(filename, tol)

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

Fs=1/(time(2)-time(1));

base_accel=base_accel-mean(base_accel);
base_disp=base_disp-mean(base_disp);
floor_1=floor_1-mean(floor_1);
floor_2=floor_2-mean(floor_2);
floor_3=floor_3-mean(floor_3);
floor_4=floor_4-mean(floor_4);

figure
plot(time,base_disp);
grid on
title('Base disp');

figure
plot(time,base_accel);
grid on
title('Base Accel - Data');
xlabel('time');
ylabel('acceleration');

figure
Y = fft(base_accel);
L=length(time);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

[pk loc]=findpeaks(P1,'MinPeakHeight',tol);
amp_array1=pk;
freq_array1=f(loc);

invertedP1=max(P1)-P1;
[pk loc]=findpeaks(invertedP1,'MinPeakHeight',tol);
amp_array_valley=max(P1)-pk;
freq_array_valley=f(loc);
counter=0;
freq_array2=0;
for i=1:length(amp_array_valley)
    if amp_array_valley(i)>=tol
        counter=counter+1;
        amp_array2(counter)=amp_array_valley(i);
        freq_array2(counter)=freq_array_valley(i);
    end
end


amp_array=amp_array1;
freq_array=freq_array1;
hold on
grid on
plot(f,P1) 
plot(freq_array1,amp_array1,'ro');
% if freq_array2~=0
%     plot(freq_array2,amp_array2,'r+');
% end
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


figure
NFFT=2^(ceil(log2(length(base_accel))));
L=length(base_accel);         
X=fftshift(fft(base_accel,NFFT));         
Px=X.*conj(X)/(NFFT*L); %Power of each freq components       
fVals=Fs*(-NFFT/2:NFFT/2-1)/NFFT;
hold on
plot(fVals,Px,'b');  
title('One Sided Power Spectral Density');
xlabel('Frequency (Hz)')
ylabel('PSD');
grid on
