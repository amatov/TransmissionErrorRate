function simustep_old(N)

%HowTo
%-----
%We consider here 3 cables: At the end of the first the both(l0), the second (l1) and the third (l2)are connected simultainusely.
%Maindatainputs:        Zc0 = caracteristic impedance of startcable (l0)
%						Zc1 = caracteristic impedance of long cable (l1)
%						Zc2 = caracteristic impedance of short cable(l2)
%						l0  = lenght in meters of startcable (l0)
%						l1  = lenght in meters of long cable (l1)
%						l2  = lenght in meters of short cable(l3)
%						Zt1 = termination impedance of long cable (l1)
%						Zt1 = termination impedance of short cable(l2)
%						N = Number of calculationpoints used for time- and distancedimension of matrix
Zc0=27;
Zc1=27;
Zc2=45;
Zt0=50;%starimp
Zt1=12.5;
Zt2=25;
l0=10;
l1=10;
l2=30;

T=1e-6; % micros, total time for the simulation
t1=0.125e-6;
dT=T/N;
t=linspace(0,T,N);
%constant values
%c0=299792458;   

%reflection factors of lines
%start=1-(generatorvoltage/((startimp/Zc0)+1));

k=Zc0/(Zc0+Zt0);
refl0=-(Zc0-Zt0)/(Zc0+Zt0);%start/entry impedance reflection
refl1=-(Zc1-Zt1)/(Zc1+Zt1);%termination long cable
refl2=-(Zc2-Zt2)/(Zc2+Zt2);%termination short cable
ramrefl0=-(Zc0-(Zc1*Zc2/(Zc1+Zc2)))/(Zc0+(Zc1*Zc2/(Zc1+Zc2)));%reflectioncoefficient from startcable back to startcable
ramrefl1=-(Zc1-(Zc0*Zc2/(Zc0+Zc2)))/(Zc1+(Zc0*Zc2/(Zc0+Zc2)));%reflectioncoefficient from long cable back to long cable
ramrefl2=-(Zc2-(Zc1*Zc0/(Zc1+Zc0)))/(Zc2+(Zc1*Zc0/(Zc1+Zc0)));%reflectioncoefficient from short cable back to short cable

c0=299792458;%300000000; % m/s
eps=2.05;%2.05;%given by manufacurer between 2 and 2.5
v=c0/eps; %propagation velocity in lines
dl=v*dT;  % the distnace at which the signal propagates for time dT 

O0=zeros([N 1]); O1=zeros([N 1]); O2=zeros([N 1]); O3=zeros([N 1]); % the input, the junction, and the two outputs

N0=round(l0/dl); N1=round(l1/dl); N2=round(l2/dl); % the number of slices for each cabel
f0=zeros([N0 1]); b0=zeros([N0 1]); % buffer for propagation for each cable (forward&backward)
f1=zeros([N1 1]); b1=zeros([N1 1]);
f2=zeros([N2 1]); b2=zeros([N2 1]);

input=2*ones([N 1]);
o0=zeros([N 1]); om=zeros([N 1]); o1=zeros([N 1]); o2=zeros([N 1]);
for i=1:N
%   if i*dT>t1
%      input(i)=0; %input(i-1)-1;
%   end
   %if input(i)<0
   %   input(i)=0;
   %end
   i0=mod(i,N0)+1; i1=mod(i,N1)+1; i2=mod(i,N2)+1; 
   j0=mod(i+1,N0)+1; j1=mod(i+1,N1)+1; j2=mod(i+1,N2)+1; 
   
   f0(i0)=k*input(i)+refl0*b0(j0);
   b0(i0)=ramrefl0*f0(j0)+(1-ramrefl1)*Zc2/(Zc2+Zc1)*b1(j1)+(1-ramrefl2)*Zc1/(Zc2+Zc1)*b2(j2);
   
   f1(i1)=ramrefl1*b1(j1)+(1-ramrefl0)*Zc2/(Zc2+Zc0)*f0(j0)+(1-ramrefl2)*Zc0/(Zc0+Zc2)*b2(j2);
   b1(i1)=refl1*f1(j1);
   
   f2(i2)=ramrefl2*b2(j2)+(1-ramrefl0)*Zc1/(Zc0+Zc1)*f0(j0)+(1-ramrefl1)*Zc0/(Zc0+Zc1)*b1(j1);
   b2(i2)=refl2*f2(j2);
   
   o0(i)=f0(i0)+b0(j0); 
   %om(i)=f0(j0)+b0(i0)+b1(j1)+b2(j2); 
   o1(i)=(1)*(f1(j1)+b1(i1)); 
   o2(i)=(1)*(f2(j2)+b2(i2));
end  

%fid = fopen('D:\matlaba\input.dat','w');
%fprintf(fid,'(Waveform\n (numDims 1)\n (size %d)\n (dim 1\n  (extent 0 1.000000000000001E-006)\n )\n (data \n [',(N));
%fprintf(fid,'%12.8f',o0);
%fprintf(fid,' ]\n )\n)');
%fclose(fid);
c=0.005e-6; %cable constant
e=exp(-t/c);
C0=7/40*conv(e,o0);
C1=4/24*conv(e,o1);
C2=35/170*conv(e,o2);
%k=1;
%while o1(k) == 0.0 ,
%   k=k+1;
%end
%O1=o1(k-1:N);

figure(1);
%plot(linspace(0,T,N),o0,'g');%input cable 0 rectangular
[a,b,c,d]=textread('D:\MATLABa\50bignewY12.5_25_2.txt','%f %f %f %f');
hold on;
plot(a,b,'r');
plot(linspace(0,T,N),C0(1:N),'b');%input cable 0 rectangular
xlabel('Time [s]');
ylabel('Amplitude [V]');
title('Step Responce at Input');
legend('Model', 'Measurement');
hold off;

figure(2);
%plot(linspace(0,T,N),o1,'g');%end cable 1 rectangular
hold on;
plot(a,c,'r');
plot(linspace(0,T,N),C1(1:N),'b');%input cable 1 rectangular
xlabel('Time [s]');
ylabel('Amplitude [V]');
title('Step Responce at Ouput 1');
legend('Model', 'Measurement');
hold off;


figure(3);
%plot(linspace(0,T,N),o2,'g');%end cable 2 rectangular
hold on;
plot(a,d,'r');
plot(linspace(0,T,N),C2(1:N),'b');%input cable 2 rectangular
xlabel('Time [s]');
ylabel('Amplitude [V]');
title('Step Responce at Ouput 2');
legend('Model', 'Measurement');
hold off;

figure(4);
plot(linspace(0,T,N),e,'b');%input cable 0 rectangular
%hold on;
%plot(linspace(0,T,N),o1,'r');%end cable 1 rectangular
%plot(linspace(0,T,N),o2,'g');%end cable 2 rectangular
%xlabel('Time [s]');
%ylabel('Amplitude [V]');
%title('Step Response');
%legend('input','ouput1','output2');
%hold off;

