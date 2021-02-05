%Design Project Part 1
%Lillian Anderson and Bradley Blaho

offbox = readtable('C:\Users\bblah\Documents\School\Fall_2019\ECE 4445\Design Project\files\offBox.dat');
onbox = readtable('C:\Users\bblah\Documents\School\Fall_2019\ECE 4445\Design Project\files\onBox.dat');

%% Small Signal Parameter Calculations
%Initialize variables
%off box
Re = 5.321;
Rmax = 109.7091;
fs = 61.266;
ws = 2*pi*fs;
f1 = 44.4937;
f2 = 81.8295;

%Small signal parameter calculations
Res = Rmax - Re;
R1  = sqrt((Re+Res)*Re);
Qms = (fs/(f2-f1))*(sqrt((Re+Res)/Re));
Qes = Qms * (Re/Res);
Qts = (Qes*Qms)/(Qes+Qms);

%on box small signal parameter calculations
%initialize variables
Vt      = 0.02832; %meters^3
Rmax_2  = 108.6832;
Res_2   = Rmax_2 - Re;    
R1_2    = sqrt((Re + Res_2) * Re);
fc      = 89.6606; %Hz
f1_2    = 71.3465; %Hz
f2_2    = 107.6429; %Hz
Qms_2   = fc / (f2_2 - f1_2) * sqrt((Re + Res_2) / Re);
Qect    = Qms_2 * (Re/Res_2);
Vas     = Vt*((fc/fs)*(Qect/Qes) - 1);

%% Driver Lumped Element Parameter Calculations
%Open box driver parameter calculations
po = 1.18; 
c = 343; %345m/s
a = 0.1; %10cm
Sd = pi*(a^2);
Cas = Vt/(po*(c^2)); 
Mas = 1/(Cas*(2*pi*fs)^2);
Ras = 1/Qms*sqrt(Mas/Cas); 
Rae = 1/Qes*sqrt(Mas/Cas);
Cms = Cas/(Sd^2); 
Mms = 1/(Cms*(2*pi*fs)^2);
Rms = Ras*Sd^2; %Rms = 1/Qms*sqrt(Mms/Cas);
Bl2 = Re/Qes*sqrt(Mms/Cms); 
Mad = Mas - (2*8*po)/(3*pi^2*a);
Mmd = Mad*Sd^2;
Ma1 = (8*po)/(3*(pi^2)*a);
Mas = Mad + 2*Ma1;

%complex number vector
s = 1i*[offbox.Frequency]*(2*pi);

%% Motional Impedance Calculations

%DONT SCALE AT ALL MOTIONAL IMPEDANCE IS GOOD
%corrections
%scale the peak
%Bl2 = Bl2*(1.03);
%Ras = (1 - .15) * Ras;

%increase motional impedance bandwidth%Ma1 = (1 - .145) * Ma1;
%Mad = (1 - .145) * Mad;
%Cas = (1 + .17) * Cas;

%calculation of the motional impedance
Zmot = (Bl2/Sd^2) ./ (((Mad).*s) + Ras + (1./(Cas.*s)) + (2*Ma1 .* s));


%% Electrical Inductance Calculation
w2 = 20000; w2 = w2*2*pi; %set and scale to rad/s
z2_mag = 52.9897;
z2_phase = 53.6516; z2_phase = z2_phase*pi/180; %set and scale to rad

w1 = 979.9467;
z1_mag = 9.6949;
z1_phase = 31.5266;

w1 = 1944.8698;
z1_mag = 13.6899;
z1_phase = 42.2054;

%w1 = 574.9993;
%z1_mag = 7.8783;
%z1_phase = 22.1946; 

%convert
w1 = w1*2*pi; %set and scale to rad/s
z1_phase = z1_phase*pi/180; %set and scale to rad

%create a complex number
z1 = z1_mag*cos(z1_phase) + 1i*z1_mag*sin(z1_phase);
z2 = z2_mag*cos(z2_phase) + 1i*z2_mag*sin(z2_phase);

%calculate the admittance of each point
y1 = 1/z1;
y2 = 1/z2;

%lossy inducantance calculations
a = real(y1)/real(y2);
b = w2/w1;

n = log(a)/log(b);

Le = cos(n*pi/2)/(w1^n * real(y1));
LE = -1/((imag(y1) + sin(n*pi/2)/w1^n*Le)*w1);

%scale inductance values
%move the inductance up and down
Le = (1 + 1.2) * Le;
LE = (1 + 1.2) * LE;

%adjust the slope
n =  (1 - 0.07) * n;

%calculate the electrical inductance
zLe = (s).^(n)*Le;
zLE = s.*LE;
zL = (((s).^(n).*Le) .* (s.*LE))./(((s).^(n) .* Le) + (s.*LE));

%% Total Electrical, Mechanical, and Acoustic Impedance
zVC = Re + zL + Zmot;

%plot everything

freq = [offbox.Frequency];

%plot magnitude of impedance
figure;
hold on;
plot(freq, abs(zVC));
set(gca, 'YScale', 'log', 'XScale', 'log');
axis([freq(1), freq(end), 0, 120]);

plot(freq, offbox.Imp_Mag);
set(gca, 'YScale', 'log', 'XScale', 'log');
axis([freq(1), freq(end), 0, 120]);
title('Voice Coil Impedance Magnitude (Measured vs. Modeled)');
xlabel('Frequency (Hz)');
ylabel('Impedance Magnitude (Ohms)');
legend({'Modeled','Measured'},'Location','southeast');
hold off;

%{
%plot angle of impedance
figure;
hold on;
zVC_angle = angle(zVC).*180./pi; %convert to degrees
plot(freq, zVC_angle);
set(gca, 'XScale', 'log');
axis([freq(1), freq(end), -90, 90]);

measured_angle = offbox.Imp_Angle;
plot(freq, measured_angle);
set(gca, 'XScale', 'log');
axis([freq(1), freq(end), -90, 90]);
title('Voice Coil Impedance Angle (Measured vs. Modeled)');
xlabel('Frequency (Hz)');
ylabel('Impedance Angle (Degrees)');
legend({'Modeled','Measured'},'Location','southeast');
%}

%% Read in excel data from multisim simulation
m_sim_data = xlsread('\\prism.nas.gatech.edu\bblaho3\ECE\ECE 4445\Design Project\3-part-circuit.xlsx');

figure;
hold on;
plot(freq, abs(zVC));
set(gca, 'YScale', 'log', 'XScale', 'log');
%axis([freq(1), freq(end), 0, 120]);
plot(m_sim_data(:,1), m_sim_data(:,2));


%% Plot the theoretical SPL
%{
%calculate the volume velocity
Bl2 = Bl2 * (1 + 6.0);
Rae = (1 - .5)*Rae;
Rat = Rae + Ras;
Ud = (Sd/Bl2) * (Rae/Rat) * (((Rat*Cas).*s))./((Mas*Cas).*(s.^2) + (Rat*Cas).*s + 1);
p = po/(2*pi) .* s .* Ud;

p = 20 * log(abs(p)/ (2e-5 * (2*pi)));


SPL = 20*log(abs(Ud) .* s / (2e-5*2*pi));

%quality factor calculation
Gs = (s / fs).^2 ./ ((s ./ fs).^2 + (1/Qts).*(s ./ fs) + 1);
p = po/(2*pi) * sqrt(Bl2)/(Sd*Re*Mas) .* Gs;
p = 20 * log(50000 * abs(p));

figure;
plot(freq, SPL);
set(gca, 'YScale', 'log', 'XScale', 'log');


%}






