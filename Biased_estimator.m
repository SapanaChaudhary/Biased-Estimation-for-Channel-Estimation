clc
clear all
close all

%% input specifications
ncb = 7e7; % nominal channel bandwidth
sf = 8e7; % sampling frequency
ts = 160e-6; % OFDM symbol duration
fd = 810; % maximum doppler freq
N=1024;

%% input data
SNR_db = 30;
SNR = 10^(SNR_db/10); % Absolute SNR
lcp = 32; % Lcp = 32
pdp = exp(-([0:lcp]/2)); % exponential PDP
P = 12; % number of pilots

%% normalization of the energy for the QPSK constelation        
temp = 0:3;           % possible symbols
temp = qammod(temp,4);  % modulated symbols
temp = abs(temp).^2;           % power of each point in the constellation
temp = mean(temp);             % average energy of the constellation
mod_norm = 1/sqrt(temp);       % normaliztion factor

%% Get the channel
% Convert variance values to linear scale
sigma= pdp;

% Find the sum of sigma and normalize 
sm = sum(sigma);
sigma = sigma/sm;

% scale the gains
for i=1:lcp
    a(i)=(sigma(i)/2)*randn(1)+1j*(sigma(i)/2)*rand(1);
end
H = (1/sqrt(N))*fft(a,1024)';
HH = H(1:108,1);
H_p = [HH(1,1) HH(2,1) HH(55,1) HH(56,1) HH(45,1) HH(46,1) HH(99,1) HH(100,1) HH(35,1) HH(36,1) HH(89,1) HH(90,1)]';
% plot(ss);
%n=10.*log10(s.^2);

%% input data generation
% symbol generation for QPSK modulation
X = randi(4,1024,1)-1;   % generation of 12x1 matrix of data symbols
% QPSK modulation
X1 = mod_norm*qammod(X,4);
X2 =(1/sqrt(N))* fft(X1,1024);
XX2 = X2(1:108,1);
Xp = [XX2(1,1) XX2(2,1) XX2(55,1) XX2(56,1) XX2(45,1) XX2(46,1) XX2(99,1) XX2(100,1) XX2(35,1) XX2(36,1) XX2(89,1) XX2(90,1)];
X_p = diag(Xp); 

% %% channel estimates from ML estimator
% Hcap_ml = zeros(5,2,1);

%% Gaussian noise covariance 
sigma2 = 1/SNR;
C_v = sigma2*eye(P);

%% Y_p is known at the output
V_p =(sigma2/2)*randn(12,1)+1j*(sigma2/2)*rand(12,1);
Y_p = X_p*H_p + V_p;

%% ML estimate of CFR at pilot locations
Cv_inv = inv(C_v);
Xp_herm = X_p';
C = C_v; %1/(X_p'*(1/C_v)*X_p); 
Hcap_ml_p = C*Xp_herm*Cv_inv*Y_p;
% Weiner smoother for Hmmse
R_hhml = xcorr2(HH,Hcap_ml_p');
R_hmlhml = xcorr2(H_p,H_p') + C;
W_opt = R_hhml*inv(R_hmlhml);
% actual CFR over the resource block : H
% requires the knowledge of channel PDP and fade characteristics
% therefore, we go uniform scattering function to design smothening filter
% vectorized optimal mmse estimate of CFR over RB
Hcap_mmse = W_opt*Hcap_ml_p;

%% Weiner smoother for using real uniform scattering function
% calculate r_rob
fdd = fd*ts; % normalized doppler freq
N = 1024;
r_rob1(1:18) = sinc(2*pi*lcp*(0:17)/N)*sinc(2*pi*fdd*(0));
r_rob2(1:18) = sinc(2*pi*lcp*(0:17)/N)*sinc(2*pi*fdd*(1));
r_rob3(1:18) = sinc(2*pi*lcp*(0:17)/N)*sinc(2*pi*fdd*(2));
r_rob4(1:18) = sinc(2*pi*lcp*(0:17)/N)*sinc(2*pi*fdd*(3));
r_rob5(1:18) = sinc(2*pi*lcp*(0:17)/N)*sinc(2*pi*fdd*(4));
r_rob6(1:18) = sinc(2*pi*lcp*(0:17)/N)*sinc(2*pi*fdd*(5));

%vectorized correlation 
r_rob = [r_rob1 r_rob2 r_rob3 r_rob4 r_rob5 r_rob6]';

% uniform Weiner smoother
W_rob = [r_rob r_rob r_rob r_rob r_rob r_rob r_rob r_rob r_rob r_rob r_rob r_rob];

%% ML estimate of CFR over the entire RB using uniform WF
Hcap_rob = W_rob*Hcap_ml_p;

%% create RB: 18 subcarriers x 6 OFDM symbols
rb = zeros(18,6,1);

%% cluster means
m11 = (Hcap_rob(1)+Hcap_rob(2)+Hcap_rob(55)+Hcap_rob(56));
m1 = 0.25*m11;
m22 = 0.25*(Hcap_rob(45)+Hcap_rob(46)+Hcap_rob(99)+Hcap_rob(100));
m2 = 0.25*m22;
m33 = 0.25*(Hcap_rob(35)+Hcap_rob(36)+Hcap_rob(89)+Hcap_rob(90));
m3 = 0.25*m33;

%% shrinkage targets
% 3 clusters -> 5 shrinkage targets 
% shrinkage target 1 
Hst_1 = (1/12)*(m11+m22+m33)*ones(12,1); 

% shrinkage target 2
Hst_2 = [m11*ones(2,1) ; m33*ones(2,1); m22*ones(2,1); m11*ones(2,1); m33*ones(2,1); m22*ones(2,1)];

%% modified JS estimate of CFR for the pilots
[temp1 temp2] = eigs(C);
C_max_eig = temp1(1,1);
p_tild = trace(C)/C_max_eig;
% second term
h2 = 1 - (p_tild - 2)/((Hcap_ml_p-Hst_2)'*inv(C)*(Hcap_ml_p-Hst_2)) ;

if (h2 >0)
    H_js = Hst_2 + h2* (Hcap_ml_p-Hst_2);
else
    H_js = Hst_2;
end

%% final CFR over the entire RB
H_js_2d =  W_rob*H_js;

%% Error
e1 = sum(conj(HH - Hcap_mmse).*(HH - Hcap_mmse)) ;
e2 = sum(conj(HH - Hcap_rob).*(HH - Hcap_rob));
e3 = sum(conj(HH - H_js_2d).*(HH - H_js_2d)) ;


