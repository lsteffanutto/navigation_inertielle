%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inertial Reference System (IRS) mechanisation in Local Tangent Plane
% (LTP).
% Simplified 2D mechanization for land vehicle.
%
% The script requires as input:
%   Reference trajectory
%   3D accelerometer measurements
%   3D gyroscope measurements
%   Initial attitude
%   Initial position
%   Initial velocity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright a ENAC, 2021
% ENAC : http://www.enac.fr/
% signav@recherche.enac.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear ; 
close all; clc; format long
addpath(genpath('./Library'));

iFigure = 1;

%% Download data
% =========================
% Choice of trajectory
bTrajectory = 1;
DataPreprocessing(bTrajectory);
load ('Reference.mat');
load ('IMU.mat');

%% Pre-process IMU measurements
% =========================
% Specific force along platform x-axis (p-frame)
% ---------------------------------------
close(figure(iFigure)); figure(iFigure); 
plot(t_IMU, f_bi_p(:,1), 'b', 'LineWidth',2);
xlabel('Time [s]'); ylabel('Specific force [m/s2]');
grid on; title('Specific force measurement along platform x-axis, f_b_/_i ^p(x)');
iFigure=iFigure+1;
% Angular rate along platform z-axis (p-frame)
% ---------------------------------------?
close(figure(iFigure)); figure(iFigure);
plot(t_IMU, w_bi_p(:,3), 'b', 'LineWidth',2);
xlabel('Time [s]'); ylabel('Specific force [m/s2]');
grid on; title('Gyrometer measurement along platform z-axis, \omega_b_/_i ^p(z)');
iFigure=iFigure+1;

% IMU data in body frame
% ---------------------------------------
R_platform2vehicule = [1 , 0 , 0; 0, -1, 0; 0, 0, -1];
f_bi_b = (R_platform2vehicule*f_bi_p')';
w_bi_b = (R_platform2vehicule*w_bi_p')';

% Correct heading rate from initial gyro bias
% ---------------------------------------

% derivee angle 
HeadingRate_ins = w_bi_b(:,3); % Heading rate % rad/s

bgyro	= mean(w_bi_b(:,3));
HeadingRate_ins = HeadingRate_ins - bgyro;

% Correct along track acceleration from initial accelerometer bias
% ---------------------------------------
g_ref = f_g_WGS(llh_Ref(1,1),llh_Ref(1,3));
aAT_ins = f_bi_b(:,1); % Along track (AT) acceleration % m/s2

bacc = mean(f_bi_b(:,1));
aAT_ins   = aAT_ins-bacc + g_ref;

f_bi_b_x = f_bi_b(:,1);
w_bi_b_z = w_bi_b(:,3);

%% Initialize output variables
% =========================
vAT_ins   = zeros(nPts_IMU,1); % Along track (AT) velocity estimate
psi_ins   = zeros(nPts_IMU,1); % IRS heading estimate
theta_ins = zeros(nPts_IMU,1); % IRS pitch estimate
vn_ins    = zeros(nPts_IMU,1); % IRS North velocity estimate in LTP
ve_ins    = zeros(nPts_IMU,1); % IRS East velocity estimate in LTP
pn_ins    = zeros(nPts_IMU,1); % IRS position estimate in LTP - North axis
pe_ins    = zeros(nPts_IMU,1); % IRS position estimate in LTP - East axis


%% Implement IRS solution
% =========================

h = waitbar(0,'IRS navigation ...');
j=1;
for i=1:nPts_IMU
    
    if (mod(i,100) == 0)  % Test GPS data availability
        waitbar(i / nPts_IMU);
        j=j+1;
    end
    
    % simplified IRS mechanisation
    % ============================

    if (i == 1) % Initialize IRS platform
        dT = Ts;
        psi_ins(i) = heading_Ref(1);
        vAT_ins(i) = 0;
        vn_ins(i)  = 0;
        ve_ins(i)  = 0;
        pn_ins(i)  = ned_Ref(1);
        pe_ins(i)  = ned_Ref(1);

    else % IRS 2D mechanization
        dT = dT + Ts;
        psi_ins(i) = psi_ins(i-1) + Ts*(HeadingRate_ins(i)+ HeadingRate_ins(i-1))/2;
        aAT_ins(i) = f_bi_b_x(i) + aAT_ins(i-1);
        vAT_ins(i) = vAT_ins(i-1) + Ts*( aAT_ins(i) +  aAT_ins(i-1) )/2 ;
        vn_ins(i)  = vn_ins(i-1)  + ( vAT_ins(i)*cos(HeadingRate_ins(i)) + vAT_ins(i-1)*cos(HeadingRate_ins(i-1)) )/2;
        ve_ins(i)  = ve_ins(i-1)  + ( vAT_ins(i)*sin(HeadingRate_ins(i)) + vAT_ins(i-1)*sin(HeadingRate_ins(i-1)) )/2;
        pn_ins(i)  = pn_ins(i-1)  + Ts*(vn_ins(i) + vn_ins(i-1))/2;
        pe_ins(i)  = pe_ins(i-1)  + Ts*(ve_ins(i) + ve_ins(i-1))/2;
        
    end
    
end
close(h);

%% Display Results
% =========================

% Define a common time vector
% ---------------------------------------
t = 0:1/100:(nPts_IMU-1)/100;
t_corr = -(nPts_IMU-1)/100:1/100:(nPts_IMU-1)/100;
t = t';
t_corr = t_corr';
% AT acceleration
% ---------------------------------------
close(figure(iFigure)); figure(iFigure); 
plot(t, aAT_ins, 'b', 'LineWidth',2);
xlabel('Time [s]'); ylabel('Acceleration [m/s?]');
grid on; title('Along Track acceleration');
iFigure=iFigure+1;

% Heading rate
% ---------------------------------------
close(figure(iFigure)); figure(iFigure); 
plot(t, HeadingRate_ins, 'b', 'LineWidth',2);
xlabel('Time [s]'); ylabel('Angular rate [rad/s]');
grid on; title('Heading rate');
iFigure=iFigure+1;

% 2D position in LTP
% ---------------------------------------
close(figure(iFigure)); figure(iFigure); hold on; grid on; title('2D Position in LTP - NED frame');
plot(ned_Ref(:,2), ned_Ref(:,1), '.r', 'MarkerSize',24);
plot(pe_ins, pn_ins, '.-b', 'LineWidth',2, 'MarkerSize',12);
xlabel('East axis [m]'); ylabel('North axis [m]');
legend('REF','Simplified IRS');
iFigure=iFigure+1;

% Position error
% ---------------------------------------
close(figure(iFigure)); figure(iFigure); 
subplot(2,1,1), plot(t_Ref,ned_Ref(:,1)-pn_ins(1:100:end), 'g', 'LineWidth',2);
grid on,
title('IRS position error: REF - IRS');
xlabel('Time [s]'); ylabel('North position error [m]');
subplot(2,1,2), plot(t_Ref,ned_Ref(:,2)-pe_ins(1:100:end), 'k', 'LineWidth',2);
grid on,
xlabel('Time [s]'); ylabel('East position error [m]');
iFigure=iFigure+1;

% North-East velocity in LTP
% ---------------------------------------
close(figure(iFigure)); figure(iFigure); 
subplot(2,1,1),plot(t_Ref, v_ned_Ref(:,1), 'r', 'LineWidth',2);hold on; grid on; 
plot(t_Ref, vn_ins(1:100:end), 'b', 'LineWidth',2);
xlabel('Time [s]'); ylabel('North velocity [m/s]'); legend('REF','Simplified IRS');
title('Vehicle velocity in LTP');
subplot(2,1,2),plot(t_Ref, v_ned_Ref(:,2), 'r', 'LineWidth',2);hold on; grid on;
plot(t_Ref, ve_ins(1:100:end), 'b', 'LineWidth',2);
xlabel('Time [s]'); ylabel('East velocity [m/s]'); 
iFigure=iFigure+1;

% Vehicle heading
% ---------------------------------------
close(figure(iFigure)); figure(iFigure); 
subplot(2,1,1),hold on; grid on; title('Vehicle Heading');
plot(t_Ref, 180/pi*unwrap(heading_Ref), 'r', 'LineWidth',2);
plot(t_Ref, 180/pi*unwrap(psi_ins(1:100:end)), 'b', 'LineWidth',2);
xlabel('Time [s]'); ylabel('[deg]'); legend('REF','Simplified IRS');
subplot(2,1,2),hold on; grid on; title('Heading error: Reference - IRS');
plot(t_Ref, 180/pi*unwrap(heading_Ref) - 180/pi*unwrap(psi_ins(1:100:end)), 'k', 'LineWidth',2);
xlabel('Time [s]'); ylabel('[deg]');
iFigure=iFigure+1;

