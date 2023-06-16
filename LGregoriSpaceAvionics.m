%% Space Avionics Final Essay 

clear all; clc;

d2r = pi/180;
r2d = 180/pi;

%% symbol rate, information rate, BW
 
% coding_rate = m * symbol_rate %baud per sec
% information_rate = R_coding * coding_rate %bit per sec, bit rate
% occupied_BW = (1+roll_fact)*symbol_rate %Hz

%% IQ recordings 8QPSK (L-S band)
ligthspeed = 3*10^8; %m/s 
frequency_min = 1.6*10^9; %GHz
frequency_max = 2.6*10^9; %GHz
BW_recordings = 4*10^6;   %MHz

lambda_min = ligthspeed/frequency_max; %wavelength
lambda_max = ligthspeed/frequency_min; %wavelength
%bit rate, coding rate, symbol rate
m=8;
r=3/4;
roll_fact = 0.2;
Rs_recordings = BW_recordings/(1+roll_fact); %MBaud
Rc_recodings = m*Rs_recordings; %MBaud
Rb_recordings = r* Rc_recodings; %Mbit per sec
%carrier to noise ratio, SNR 
system_temp = 290; %K, noise system temperature
Boltz_cost = -228.6; %dBW/KelvinHz
atm_loss = 0.49; %dB, losses due to atmosphere  
pointing_loss = 0.55; %dB, losses due to the receiver/GS position
range = 600*10^3; % km, LEO orbit for microsatellite missions (300km-1000km)

%free space path loss 
FSP_min = 20*log10(4*pi*range/lambda_max); %dB
FSP_max = 20*log10(4*pi*range/lambda_min); %dB

%atomospheric losses
losses_min = FSP_min + atm_loss + pointing_loss; %dB
losses_max = FSP_max + atm_loss + pointing_loss; %dB

%transmitter antenna gain 
antenna_diameter = 5.1324; %meters, diameter
sat_area = pi*(antenna_diameter/2)^2;
eta_area = 0.7;
Ae_tx = eta_area * sat_area;
Gsat_min = (4*pi*Ae_tx/(lambda_max^2));
Gsat_max = (4*pi*Ae_tx/(lambda_min^2));
Gsat_mindB = 10*log10(Gsat_min);
Gsat_maxdB = 10*log10(Gsat_max);

%receiver antenna gain 
semimajor_axis = 0.77/2; %m
semiminor_axis = 0.72/2; %m
receiver_area = pi*semiminor_axis*semimajor_axis; %m2
eta_area = 0.7;
Ae_rx = eta_area * receiver_area;
Greceiver_min = (4*pi*Ae_rx/(lambda_max^2));
Greceiver_max = (4*pi*Ae_rx/(lambda_min^2));
Greceiver_mindB = 10*log10(Greceiver_min);
Greceiver_maxdB = 10*log10(Greceiver_max);

%eirp sat
Pt = 1; %watt, power transmitted to the antenna
%Gt = 50; %dB, standard value, antenna satellite gain 
eirp_sat_min = Pt*Gsat_mindB; 
eirp_sat_max = Pt*Gsat_maxdB;

SNR_min = eirp_sat_min - losses_min + Greceiver_mindB - Boltz_cost - 10*log10(system_temp*BW_recordings); %dB;%carrier to noise ratio 
SNR_max = eirp_sat_max - losses_max + Greceiver_maxdB - Boltz_cost - 10*log10(system_temp*BW_recordings); %dB;%carrier to noise ratio 

%if limited in bandwidth Shannon's Theorem
% bandwidth-limited, SNR > 0

Channel_capacity_min = 0.332*BW_recordings*SNR_min;
Channel_capacity_max = 0.332*BW_recordings*SNR_max;

%% TT&C , 8PSK (L-S band)
Rb_TTC = 38.4*10^3 ; %kbit per sec 
frequency_TTC = 2.0*10^9; %GHz
lambda_TTC = ligthspeed/frequency_TTC; %wavelength
%choose the coding parameters to derive the occupied bandwidth
m_TTC = 8 ;%number of bits per symbol
r_TTC = 3/4 ;%FEC rate, code rate
Rc_TTC = Rb_TTC/r_TTC;  %coding rate
roll_fact_TTC = 0.2; %roll factor for coding
Rs_TTC = Rc_TTC/m_TTC; %symbol rate
BW_TTC = (1+roll_fact_TTC)*Rs_TTC; %occupied bandwidth

%free space path loss 
FSP_TTC = 20*log10(4*pi*range/lambda_TTC); %dB
%atomospheric losses
losses_TTC = FSP_TTC + atm_loss + pointing_loss; %dB
%transmitter antenna gain 
Gsat_TTC = (4*pi*Ae_tx/(lambda_TTC^2));
Gsat_TTCdB = 10*log10(Gsat_TTC);
%receiver gain
Greceiver_TTC= (4*pi*Ae_rx/(lambda_TTC^2));
Greceiver_TTCdB = 10*log10(Greceiver_TTC);
%eirp sat
Pt = 1; %watt, power transmitted to the antenna
eirp_sat_TTC = Pt*Gsat_TTCdB; 
%snr ttc 
SNR_TTC = eirp_sat_TTC - losses_TTC + Greceiver_TTCdB - Boltz_cost - 10*log10(system_temp*BW_TTC); %dB;%carrier to noise ratio 
%channel capacity TTC
Channel_capacity_TTC = 0.332*BW_TTC*SNR_TTC;

%BER_TTC = 1*10^(-10); %bit error rate CONTROLLARE VALORE!!!!

%% PAYLOAD DATA 8PSK (X-Band)
Rb_payload = 20*10^6; %Mbit per sec
frequency_payload = 10.0*10^9; %GHz
lambda_payload = ligthspeed/frequency_payload; %wavelength
m_payload = 8 ;%number of bits per symbol
r_payload = 3/4 ;%FEC rate, code rate
Rc_payload = Rb_payload/r_payload ;%coding rate
roll_fact_payload = 0.2; %roll factor for coding
Rs_payload = Rc_payload/m_payload; %symbol rate
BW_payload = (1+roll_fact_payload)*Rs_payload; %occupied bandwidth

%free space path loss 
FSP_payload = 20*log10(4*pi*range/lambda_payload); %dB
%atomospheric losses
losses_payload = FSP_payload + atm_loss + pointing_loss; %dB
%transmitter antenna gain 
Gsat_payload = (4*pi*Ae_tx/(lambda_payload^2));
Gsat_payloaddB = 10*log10(Gsat_payload);
%receiver gain
Greceiver_payload= (4*pi*Ae_rx/(lambda_payload^2));
Greceiver_payloaddB = 10*log10(Greceiver_payload);
%eirp sat
Pt = 1; %watt, power transmitted to the antenna
eirp_sat_payload = Pt*Gsat_payloaddB; 
%snr ttc 
SNR_payload = eirp_sat_payload - losses_payload + Greceiver_payloaddB - Boltz_cost - 10*log10(system_temp*BW_payload); %dB;%carrier to noise ratio 
%channel capacity TTC
Channel_capacity_payload = 0.332*BW_payload*SNR_payload;
%BER_payload = 1*10^(-10); %bit error rate CONTROLLARE VALORE!!!!

%% scrub rate 
mu_earth= 3.986004418 *10^14; %gravitational earth parameter m3/s2 since m_sat<<M_earth its negligible in the total calculation. 
r=(600+6378)*10^3; %km
mean_motion = sqrt(mu_earth/r^3);
orbital_period = 2*pi/mean_motion;
dEdt = 30/(orbital_period*8*10^6) ; %error rate in a orbit converted in second, 30bit/MB/orbit
used_bits = 1;
total_bits = 10;
mit_windows = 5000;
fanout = 10; %between (1,mit_windows)
dCdt = ((dEdt)^2)*((used_bits/total_bits)^2)*(2/3)*(1/mit_windows)*fanout; %scrub rate, scrub/second

%scrub rate as 10 times the given error rate if it is based on SEU rate
scrub_rate = 10*dEdt;

%% PART3 TASKS
B_bound1 = 5*(2^(1/5)-1)
task1 = 1/8;
task2 = 1/8; 
task3 = 2/7;
task4 = 1/12;
task5 = 2/5;
cpu_Utilisation1 = (task1 + task2 + task3 + task4 + task5)

%7 tasks in 12 period of time
newtask1 = 1/8;
newtask2 = 1/8; 
newtask3 = 2/10;
newtask4 = 1/12;
newtask5 = 2/10;
cpu_Utilisation2 = (newtask1 + newtask2 + newtask3 + newtask4 + newtask5)

cpu_timereduced = cpu_Utilisation2 /0.8