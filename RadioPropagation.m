clc;
close all;
clear all;
%% task 1

d = 10:0.1:10000; %distance d from 10m to 10km
Rd = sqrt(d.^2 + 35.7^2);
f = 9*1e8; %for frequency of 900MHz
%f = 18*1e8; %for frequency of 1800MHz
%f = 35*1e8; %for frequency of 3500MHz
Er = 16.3;

Beta = asin((37+1.3)./d);
k = 2*pi*f/(3*1e8); % definition of k, wavenumber in free space 
E0 = (exp(-1i*k.*Rd))/Rd;
landa = (3*1e8)/f;
Rv = (-Er*sin(Beta)+sqrt(Er-cos(Beta).^2))/(Er*sin(Beta)+sqrt(Er-cos(Beta).^2));
%Reflection coefficient for vertically polarized wave
Rh = (sin(Beta)-sqrt(Er-cos(Beta).^2))/(sin(Beta)+sqrt(Er-cos(Beta).^2));
% Reflection coefficient for horizontally polarized wave

R1 = 37 ./ sin(Beta);
R2 = 1.3 ./ sin(Beta);
Ri = R1 + R2;
delta = Ri - Rd;
%Evh = ((norm(E0))^2)*4*(sin(2*pi*37*1.3/(landa*d)))^2; % As hT and hR are very smaller than d
Ev = abs(E0*(1+Rv*exp(-1i*k*delta))); % Electric field strength for vertcally polarized wave
Eh = abs(E0*(1+Rh*exp(-1i*k*delta))); % Electric field strength for horizontally polarized wave
db = (4*37*1.3)/landa; % The breakpoint distance
disp(['break point distance is ' num2str(db) 'm']);

figure;
semilogx(d,10*log10(Ev),'r', d, 10*log10(Eh), 'b','LineWidth',0.5);
legend('Ev', 'Eh');
xlabel('Distance')
ylabel('Electric field')
grid on;
hold on;

%% task 2
% Break point distance and E field strength for three different frequencies were displayed in
% previous part

%% task 3
f = 9*1e8  

%for frequency of 900MHz
d = 10:0.1:10000; %distance d from 10m to 10km

L_2ray = 40*log10(d) - 20*log10(37) - 20*log10(1.3); % The path loss for 2ray reflection model(dB)

L_fs = 32.4 + 20* log10(d/1e3) + 20*log10(f/1e6); % The free space path loss for 900MHz

ahms = ((1.1*log10(f/1e6)-0.7)*1.3) - (1.56*log10(f/1e6)-0.8);
L_cost231 = 69.55 + 26.16*log10(f/1e6) - 13.82*log10(37) - ahms + (44.9-6.55*log10(37))*log10(d/1e3) + 20;
% Path loss for COST 231 model

figure;
subplot(3,1,1);
semilogx(d,L_2ray,'k');
legend('Path loss for 2Ray model');
grid on;
subplot(3,1,2);
semilogx(d,L_fs,'r');
legend('Free space path loss');
grid on;
subplot(3,1,3);
semilogx(d,L_cost231,'b');
legend('Path loss for COST 231 model');
xlabel('Distance')
ylabel('Pathloss')
title('Three different propagation models with f=900MHz')
grid on;
hold on;

f = 18*1e8 %for frequency of 1800MHz
L_2ray = 40*log10(d) - 20*log10(37) - 20*log10(1.3); % The path loss for 2ray reflection model(dB)

L_fs = 32.4 + 20* log10(d/1e3) + 20*log10(f/1e6); % The free space path loss for 1800MHz

ahms = ((1.1*log10(f/1e6)-0.7)*1.3) - (1.56*log10(f/1e6)-0.8)
L_cost231 = 46.3 + 33.9*log10(f/1e6) - 13.82*log10(37) - ahms + (44.9-6.55*log10(37))*log10(d/1e3) + 20;
% Path loss for COST 231 model

figure;
subplot(3,1,1);
semilogx(d,L_2ray,'k');
legend('Path loss for 2Ray model');
grid on;
subplot(3,1,2);
semilogx(d,L_fs,'r');
legend('Free space path loss');
grid on;
subplot(3,1,3);
semilogx(d,L_cost231,'b');
legend('Path loss for COST 231 model');
xlabel('Distance')
ylabel('Pathloss')
title('Three different propagation models with f=1800MHz')
grid on;
hold on;

f = 35*1e8 %for frequency of 3500MHz
L_2ray = 40*log10(d) - 20*log10(37) - 20*log10(1.3); % The path loss for 2ray reflection model(dB)

L_fs = 32.4 + 20* log10(d/1e3) + 20*log10(f/1e6); % The free space path loss for 3500MHz

ahms = ((1.1*log10(f/1e6)-0.7)*1.3) - (1.56*log10(f/1e6)-0.8)
L_cost231 = 46.3 + 33.9*log10(f/1e6) - 13.82*log10(37) - ahms + (44.9-6.55*log10(37))*log10(d/1e3) + 20;
% Path loss for COST 231 model

figure;
subplot(3,1,1);
semilogx(d,L_2ray,'k');
legend('Path loss for 2Ray model');
grid on;
subplot(3,1,2);
semilogx(d,L_fs,'r');
legend('Free space path loss');
grid on;
subplot(3,1,3);
semilogx(d,L_cost231,'b');
legend('Path loss for COST 231 model');
xlabel('Distance')
ylabel('Pathloss')
title('Three different propagation models with f=3500MHz')
grid on;
hold on;







