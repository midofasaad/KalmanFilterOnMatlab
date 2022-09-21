clear all
close all;

global H R Q A;
%% Measurement Emulation 
a0=1;
b0=20;
% [3 1 4 2 1].'
% [0 0.1 0.5 0.6 1].'

%% Load parameters & Current from Mat file
% Idt=load('current_parameters.mat').Idt;
% a_mat=load('current_parameters.mat').a;
% b_mat=load('current_parameters.mat').b;
% c_mat=load('current_parameters.mat').c;

%% Voltage Signal
fs=1000;
delta_t=1/fs;
size_bf=1000;
tmax=size_bf*delta_t;
freq=100;
omega=2*pi*freq;
Amp_I=20 ;
R_on=3;
Amp_V=R_on*Amp_I; 
phase_shift=20;


t=linspace(0,tmax,size_bf);
x_ref=Amp_I*sin(omega*t+phase_shift);
x_ref_2=Amp_V*sin(omega*t+phase_shift);




%% Voltage Measurment Model parameters

omega_k=omega/fs;
gain=1;
sigma=20;
variance=sigma^2;
offset=20;
noise = offset+ sigma*rand(1,length(x_ref));

y=x_ref+noise;
y_2=x_ref_2+noise; 

y_kal=[y;y_2];
% Yw_test=fft(y)



%% Gaussian Recursive estimation parameters
A_estimate=1;
A_estimates=[];
phase_estimate=2;
phase_estimates=[];
gaus_estimates=[];
c=1;
lamda=0.95;


%% Kalman Filter intial factors and Matrices

%Prediction Matrices & factors
A_c=[0,1,0,0,0,0;
    -omega^2,0,0,0,0,0;
    0,0,0,1,0,0;
    0,0,-omega^2,0,0,0;
    0,0,0,0,0,0;
    0,0,0,0,0,0];

X_kal =[1;1;1;1;1;1];
P =[10,0,0,0,0,0;
    0,10,0,0,0,0;
    0,0,10,0,0,0;
     0,0,0,10,0,0;
     0,0,0,0,10,0;
     0,0,0,0,0,10];
A =[cos(omega_k),sin(omega_k)/omega,0,0;
    -sin(omega_k)*omega,cos(omega_k),0,0;
    0,0,1,0;
    0,0,0,1];

Q = [0,0,0,0,0,0;
    0,0,0,0,0,0;
    0,0,0,0,0,0;
     0,0,0,0,0,0;
     0,0,0,0,0,0;
     0,0,0,0,0,0];

B = [0;0;0;0;0;0];
U = 0;

% Correction Matrices

H = [0,1,1,0,0,0;
    0,0,0,1,0,0];
R = [variance,0 ;
    0,variance]; 
% c_system=ss(A_c,B,H,0);
% d_system=c2d(c_system,delta_t);
% Ad=d_system.A;
%Storage
kalman_estimates=[];
variance=[];

%% Particle Filter Parameters and intializations
a=-20
b=20
N=1000
particles=a + (b-a)*rand(2,N); 
weights=ones(2,N)/N ;
z=delta_t*omega;
r=1/((z)^2+1);
x_pf=[];
x_pf_variance=[];

%% Berto-yoshida intialization parameters
online_measurements=[];
yoshida_estimate=[]



%% Offset and Variance estimation intialliazions
variance_estimates=[];
sigma_estimates=[];
offset_estimate=0;
offset_estimates=[];
diff_square=0;
variance_estimate=0;


%% Estimation Loop
for k=1:length(y)
 % Offset estimation and compensation
     offset_estimate=(y(k)+(k-1)*offset_estimate)/k;
     offset_estimates=[offset_estimates, offset_estimate];
     y_corr=y(k)-offset_estimate;
 %Bertocco-yoshida estimation 
     online_measurements=[online_measurements,y(k)];
 if k>2
     Yw=fft(online_measurements);
     [Yabs, id] = max(abs(Yw)); %why half of the spectrum?
     
     k_before = id; 
     k_peak=id; 
     k_after = id+1; 

     w_k = 2*pi/k; 
     
     wkm1= (k_before-1)*w_k;
     wk = (k_peak -1)*w_k; 
     wkp1= (k_after-1)*w_k;
     
     if id==1
     wkm1=0 ;
     R_yoshida = (-Yw(k_peak))/(Yw(k_peak)-Yw(k_after)); 
     else 
     R_yoshida = (Yw(k_before)-Yw(k_peak))/(Yw(k_peak)-Yw(k_after)); 
     end 
     
     r_yoshida = ( -exp(-1i*wk)+exp(-1i*wkm1) )/( -exp(-1i*wkp1)+exp(-1i*wk) ); 
     lambda = exp(1i*wk)*(r_yoshida-R_yoshida)/( r_yoshida*exp(-1i*2*pi/k_peak)-R_yoshida*exp(1i*2*pi/k_peak) ); 
     
     c_yoshida = (1-lambda^k)/(1-lambda*exp(-1i*wk)); 
     c_yoshida = Yw(k_peak)/c; 
     A_yoshida = abs(c); 
     p_yoshida = angle(c);
     
     yoshida_estimate=A_yoshida*sin((omega/fs)*k+p_yoshida);
 end   
        
 
     
 %Gaussian Estimation
     c=lamda*c+0.5;
     e=y_corr-A_estimate*sin((omega/fs)*k+phase_estimate);
     phase_estimate_prev=phase_estimate;
     phase_estimate=phase_estimate+cos((omega/fs)*k+phase_estimate)*(e/(A_estimate*c));
     A_estimate=A_estimate+sin((omega/fs)*k+phase_estimate_prev)*(e/c);
     gaus_estimate=A_estimate*sin((omega/fs)*k+phase_estimate);
     
     phase_estimates=[phase_estimates,phase_estimate];
     A_estimates=[A_estimates,A_estimate];
     gaus_estimates=[gaus_estimates,gaus_estimate];
    
%Linear Kalman Filter
     %Prediction step
     F=[cos(omega_k),sin(omega_k)/omega,0,0,0,0;
    -sin(omega_k)*omega,cos(omega_k),0,0,0,0;
    0,0,cos(omega_k),sin(omega_k)/omega,0,0;
    0,0,-sin(omega_k)*omega,cos(omega_k),0,0;
    0,0,0,0,1,0;
    delta_t*(2*X_kal(2)*X_kal(3)/X_kal(1)^3 -X_kal(4)/X_kal(1)),-delta_t*X_kal(3)/X_kal(1)^2,delta_t*(-X_kal(2))/X_kal(1)^2,delta_t/X_kal(1),0,1];

     X_kal =[cos(omega_k)*X_kal(1)+(sin(omega_k)/omega)*X_kal(2) ;
             -sin(omega_k)*omega*X_kal(1)+cos(omega_k)*X_kal(2);
             cos(omega_k)*X_kal(3)+(sin(omega_k)/omega)*X_kal(4);
             -sin(omega_k)*omega*X_kal(3)+cos(omega_k)*X_kal(4);
             X_kal(5);
             X_kal(6)+delta_t*(X_kal(1)*X_kal(4)-X_kal(2)*X_kal(3))/X_kal(1)];


     P = F*P*F' + Q;
     % Update H matrix 
      Gm=[1,0,0,0,1,0;
         0,0,1,0,X_kal(6),X_kal(6)];

     ym = [X_kal(1)+X_kal(5); X_kal(3)+X_kal(6)*X_kal(5)];
     Pm = Gm*P*Gm' +R;
     K = P*Gm'*inv(Pm);       
     X_kal = X_kal + K*(y_kal(:,k)-ym);
     P = P - K*Gm*P;
     

     
     kalman_estimates=[kalman_estimates, X_kal];
     variance=[variance, P(1,1)];
%Particle Filter
    % Prediction Step
%     particles=A*particles;
     
     
    % Correction Step
    
%variance estimation 

     variance_estimate=(k-1)/(k)*variance_estimate+(1/k)*(y(k)-offset-x_ref(k))^2;
     sigma_estimate=sqrt(variance_estimate);
     
     variance_estimates=[variance_estimates, variance_estimate];
     sigma_estimates=[sigma_estimates, sigma_estimate];
    
end 

% P =[1000,0,0;
%     0,1000,0;
%     0,0,1000];



gaus_RMSE=RMSE(gaus_estimates,x_ref);
kalman_RMSE=RMSE(kalman_estimates(2,:),x_ref);


%% Plotting area
% labels={'Original signal','PF_Estimate','KF_Estimate','Gaus_Estimate','Hilbert_Estimate'}
labels={'Original signal','Measurement','Gaus-Estimate','KF-Estimate'}

figure()
plot(t,x_ref,'red');
hold on 
plot(t,y,'blue');
plot(t,gaus_estimates,'green');
plot(t,kalman_estimates(1,:),'black');
grid('on')
legend(labels)

figure()
plot(t,gaus_RMSE,'green');
hold on 
plot(t,kalman_RMSE,'black');
grid('on')
legend('gaus-RMSE','KF-RMSE');

figure()
plot(t,R_on);
hold on 
plot(t,kalman_estimates(4,:))


