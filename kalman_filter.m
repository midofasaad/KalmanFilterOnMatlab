clear all
close all;

%% Load Base Parameters for the Lookup table
load('Lookup.mat')
load('current_parameters.mat')
load('parameters')


% Testing Lookup tables
current=200; 
R_on_meas=2; %mohm
[Ron_k, Tj_k]= determine_fitting_points(Idt,Tj,Rds_on,current) 
[Tj,parameters]=determine_Tj(R_on_meas,Tj_k,Ron_k)
Tj_3= determine_Tj_3(parameters_a,parameters_b,parameters_c,current,R_on_meas)


%% Measurement Emulation 
I_sig=1;
I0_off=20;
Ids_gain=10;
I_var=I_sig^2;

V_sig=1;
V0_off=0.2;
Vds_gain=1;
V_var=V_sig^2;

Qi = [5,0,0;
     0,5,0;
     0,0,10];
 
% Q_2 = [0,0,0,0;
%     0,0,0,0;
%     0,0,0,0;
%     0,0,0,0];
Q_2 = [5,0,0,0,0,0;
    0,5,0,0,0,0;
    0,0,5,0,0,0;
    0,0,0,5,0,0;
    0,0,0,0,10,0;
    0,0,0,0,0,10];

 
Qv = [5,0,0;
     0,5,0;
     0,0,10];
 
 

%% Voltage Signal
fs=1000;
delta_t=1/fs;
freq=20;
omega=2*pi*freq;
omega_k=omega/fs;
gain=1;

%% Kalman Filter intial factors and Matrices

%Prediction Matrices & factors
A_c=[0,1,0;
    -omega^2,0,0;
    0,0,0];
A_c_2=[0,1,0,0;
    -omega^2,0,0,0;
    0,0,0,0;
    0,0,0,0];

X_kal =[0;0;0];
X_kal_2=[0;0;0;0;1;2];

P =[10,0,0;
    0,10,0;
    0,0,10];
P_2 =[300,0,0,0;
    0,300,0,0;
    0,0,300,0;
     0,0,0,300];

P_3=[100,0,0,0,0,0;
    0,100,0,0,0,0;
    0,0,100,0,0,0;
    0,0,0,100,0,0;
    0,0,0,0,100,0;
    0,0,0,0,0,100];

A =[cos(omega_k),sin(omega_k)/omega,0,0;
    -sin(omega_k)*omega,cos(omega_k),0,0;
    0,0,1,0;
    0,0,0,1];
B = [0;0;0];
B_2=[0;0;0;1];

Bv = [1;0;0];
U = 0;

% Correction Matrices

H = [1,0,1];
H_2 = [1,0,0,0,1,0;
    0,0,1,0,0,1];
Rv = V_var ;
Ri = I_var ;
R_2=[I_var, 0;
    0,V_var];   

% c_system=ss(A_c,B,H,0);
% c_system_2=ss(A_c_2,B_2,H_2,0);
% d_system=c2d(c_system,delta_t);
% d_system_2=c2d(c_system_2,delta_t);
% Ad=d_system.A;
% Ad_2=d_system_2.A;

%% Gaussian Recursive estimation parameters
A_estimate=1;
A_estimates=[];
phase_estimate=2;
phase_estimates=[];
gaus_estimates=[];
c=1;
lamda=0.95;

