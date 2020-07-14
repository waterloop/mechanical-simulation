% SLIM Design Using a Matlab Program
clear all;
clc;

%Assign necessary constants and parameters
mu0 = 4*pi*10^-7;  %permeability of air
rhow = 1.72*10^-8; %volume resistivity of copper
rhor = 2.82*10^-8; %volume resistivity of aluminum
btmax = 1.6;       %max tooth flux density
bymax = 1.3;       %max yoke flux density
J1 = 4e6;          %rotor current density

%Assign desired values for certain variables
d = 0.00635;        %aluminum thickness/rail conductor thickness
m = 3;              %number of phases
Vline = 160;        %line voltage
V1 = Vline/sqrt(3); %phase voltage
f = 120;            %frequency 
p = 4;              %number of poles
q1 = 1;             %number of slots per pole per phase
Srated = 0.1;       %rated slip/ slipping percent
Ws = 0.005;         %width of the stator/ core width

%Data from the PCP design procedure
Fsprime = 30;       %target thrust
Vcrated= 40;        %rated rotor velocity

Vs= Vcrated/(1 - Srated);   %synchronous velocity
tau = Vs/(2*f);             %pole pitch
lambda = tau/(m*q1);        %slot pitch
Ls = p*tau;                 %length of one stator unit

for i = 1:30
 N1 = p*q1*i;               %number of turns per phase
 ncos0 = 0.117930759127764; 
 ncos1(i) = 1;
 while abs(ncos0 - ncos1(i))>0.0001

     I1prime = (Fsprime*Vcrated)/(m*V1*ncos0);  %estimated rated stator RMS current
     Aw = I1prime/J1;                           %differs
     As = (10*i*Aw)/7;                          %differs
     ws = lambda/2;                             
     wt = ws;
     hs = As/ws;
     gm = 0.003;                                %mech gap
     go = gm + d;                               %magnetic air gap
     gamma = (4/pi)*(((ws/(2*go))*atan(ws/(2*go))) - log(sqrt(1 + ((ws/(2*go))^2))));
     kc = lambda/(lambda - gamma*go);
     ge = kc*go;
     kw = sin(pi/(2*m))/(q1*sin(pi/(2*m*q1)));
     G = 2*mu0*f*tau^2/(pi*(rhor/d)*ge);

     %corrections following original spreadsheet calculations
     a = ((hs/(3*ws))*(1+(3/p)))+((5*ge/ws)/(5+(4*go/ws)));
     R1(i) = 2*rhow*(pi+tau)*J1*N1/I1prime;
     X1(i) = 8*mu0*pi*f*((a*pi/q1)+(0.6*tau))*(N1^2)/p;
     Xm(i) = (24*mu0*pi*f*kw*N1^2*tau)/(pi^2*p*ge);
     R2(i) = Xm(i)/G;

     Z(i)=R1(i)+j*X1(i)+((j*R2(i)*Xm(i))/Srated)/((R2(i)/Srated) + j*Xm(i));

     I1(i) = V1/abs(Z(i)); %calculated 
     I2(i) = j*I1(i)*Xm(i)/(R2(i)/Srated+j*Xm(i));
     Im(i) = I1(i) - I2(i);

     %Actual TLIM Thrust
     Fs(i) = (m*abs(I1(i))^2*R2(i))/(((1/(Srated*G)^2)+ 1)*Vs*Srated);
     diff(i) = Fs(i) - Fsprime;
     dmin = min(abs(diff));
     Pout = Fs*Vcrated;
     Pin=Pout+m*abs(I2(i))^2*R2(i)+m*abs(I1(i))^2*R1(i);
     eta = Pout/Pin;
     PF = cos(angle(Z(i)));
     ncos1(i)=eta*PF; %ncosPhi
     ncos0=(ncos0+ncos1(i))/2; %ncosphiNew=(ncosPhiOld+ncosPhiCalculated)/2
 end;
end;
k = 1; 

while dmin~=abs(diff(k)) %checks inequality > 1 if true, 0 if false
 k = k + 1;
end;

%Step 18 question, if not we increase Nc by 1 and repeat from step 7
%works because it satisifies the condition where dmin is approximately abs
%of diff(k)

Nc = k;             %number of turns per slot, Nc
N1 = p*q1*Nc;       %number of turns per phase, N1
Fs = Fs(k);         %Thrust Force
I1 = I1(k);         %Input Current
ncos1 = ncos1(k);   %HELP > loops are confusing

%Appendix A Wire gauge Table for Copper Wire
%AWG | Diameter in mm
%{A=[3 5.8267;
%   4 5.1892;
%   5 4.6202;
%   6 4.1148;
%   7 3.6652;
%   8 3.2639;
%   9 2.90576;
%   10 2.58826]; 

A = [14 1.62814;
    15 1.45034;
    16 1.29032;
    17 1.15062;
    18 1.02362;
    19 0.91186;
    20 0.8128;
    21 0.7239]; 
gauge=0;

while (gauge<8) %8 selection of wires in the [A], need to fix
 gauge=gauge+1;
 Np = 0;
% r=0;
 wt = 1;
 wtmin= 0;
% g=0;r=0;
 while (wt-wtmin)>0.0152 %what is this condition?
%    r=r+1;
%    g=g+1;
     Dw=A(gauge,2);
     Np = Np + 1;
     ws = (Dw*10^-3*Np); %not taking into account the slot insulation factor 
     wt = lambda - ws;
     Aw = Np*pi/4*Dw^2*1e-6;
     As = (10*Nc*Aw)/7;
     hs = As/ws;
     gm = 0.003;
     go = gm + d;
     gamma=(4/pi)*(((ws/(2*go))*atan(ws/(2*go)))-log(sqrt(1 + ((ws/(2*go))^2))));
     kc = lambda/(lambda - gamma*go);
     ge = kc*go;
     G = 2*mu0*f*tau^2/(pi*(rhor/d)*ge);
     kw=sin(pi/(2*m))/(q1*sin(pi/(2*m*q1)));
     
     a = ((hs/(3*ws))*(1+(3/p)))+((5*ge/ws)/(5+(4*go/ws)));
     R1 = 2*rhow*(pi+tau)*J1*N1/I1prime;
     X1 = 8*mu0*pi*f*((a*pi/q1)+(0.6*tau))*(N1^2)/p;
     Xm = (24*mu0*pi*f*kw*N1^2*tau)/(pi^2*p*ge);
     R2 = Xm/G;

     Z=R1+j*X1+(R2/Srated*j*Xm)/(R2/Srated+j*Xm);
     I1 = V1/abs(Z);
     I2 = j*I1*Xm/(R2/Srated+j*Xm);
     Im=I1-I2;
     wtmin=2*sqrt(2)*m*kw*N1*abs(Im)*mu0*lambda/(pi*p*ge*btmax);
     J1actual=I1/Aw;
 end;
 hy=4*sqrt(2)*m*kw*N1*abs(Im)*mu0*Ls/(pi*pi*p*p*ge*bymax);
 para_wires(gauge)=Np;
 slot_width(gauge)=ws; 
 tooth_width(gauge)=wt;
 min_toothwidth(gauge)=wtmin;
 height_slot(gauge)=hs;
 Area_wire(gauge)=Aw;
 Area_slot(gauge)=As;
 Num_c(gauge)=Nc;
 Num_1(gauge)=N1;
 Sta_I(gauge)=I1;
 gap_e(gauge)=ge;
 current_den(gauge) = abs(I1)/Aw;
 height_yoke(gauge)=4*sqrt(2)*m*kw*N1*(Im)*mu0*Ls/(pi*pi*p*p*ge*bymax);
 final_thrust(gauge)=(m*abs(I1)^2*R2)/(((1/(Srated*G)^2)+1)*Vs*Srated);
 output(gauge)=final_thrust(gauge)*Vcrated;
 input(gauge)=output(gauge)+m*abs(I2)^2*R2+m*abs(I1)^2*R1;
 efficiency(gauge)= output(gauge)/input(gauge);
 difference(gauge)=final_thrust(gauge)-Fsprime;
 diffmin(gauge) = min(abs(difference));
end;

kk = min(diffmin);
jj=1;

while kk~=abs(diffmin(jj))
 jj = jj + 1;
end;
best_wiregauge=A(jj,1)
%$$$ To Generate the Characteristic curves $$$
vel_sta= Vs; 
e=1;

for slip=0.000001:0.01:1
    
 vel_rot(e)=vel_sta*(1-slip);
 if abs(Vcrated - vel_rot(e))/Vcrated < 0.01
     n_Vr = e; %index for where v = Vr
 end
 impz(e) = R1+j*X1+(R2/slip*j*Xm)/(R2/slip+j*Xm);
 i1(e) = V1/abs(impz(e));
 i2(e) = j*i1(e)*Xm/(R2/slip+j*Xm);
 im(e) = i1(e)-i2(e);
 Force(e)=(m*(abs(i1(e)))^2*R2)/(((1/(slip*G)^2)+1)*vel_sta*slip);
 out_pow(e) = Force(e)*vel_rot(e);
 in_pow(e)=out_pow(e)+m*abs(i2(e))^2*R2+m*abs(i1(e))^2*R1;
 eff(e) = out_pow(e)/in_pow(e);
 e=e+1;
end;

figure(1);
plot(vel_rot,Force,'green');
hold on;
grid on;
grid minor
ylabel('Target Thrust, Fs (N)')
xlabel('Rotor Velocity, Vr (m/s)')
plot([Vcrated Vcrated],[0,Fs])
hold on;
plot([0 Vcrated],[Fs Fs]);
hold on;
title(['Force vs. Velocity'])
legend('Actual Force','Target Velocity','Target Force')

figure(2);
plot(vel_rot,eff*100,'green');
hold on;
plot([Vcrated Vcrated],[0 eta*100]);
hold on;
plot([0 Vcrated],[eta*100,eta*100]);
hold on;
grid on
grid minor
ylabel('Efficiency (%)')
xlabel('Rotor Velocity, Vr (m/s)')
title(['Efficiency vs. Velocity'])
legend('Actual Efficiency','Target Velocity','Ideal Efficiency')