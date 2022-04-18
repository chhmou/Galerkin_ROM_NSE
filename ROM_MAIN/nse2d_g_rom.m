%{
This is a main routine for the Galerkin ROM (G-ROM)
for 2D Navier-Stokes equation

%}

clc;clear;close all

%%
format long e

global nu Re
global snapIndex
global T dt numTimeSteps

addpath ../ROM_FUN/
addpath ../FEM_FUN/


%%
%--- load DNS data and ROM data (basis and matrices)
% for simplicity, we use the ROM data generated with 1500 snapshots for
% Re=1000.

% DNS data
dns_data_dir ='/Users/cmou/Desktop/Flow_past_a_cylinder/ROM_CODE/DATA';
dns_data_name = 'nse_dns_t23_cylmesh35dt002Re_1eN3_T10.mat';
dns_data= fullfile(dns_data_dir,dns_data_name);

load(dns_data)

% POD data
load('nse_pod_re1000_snap1500.mat')

% ROM Matrices, e.g., Mr, Sr, Tr
load('nse_rom_mat_re1000_snap1500.mat')

Re = 0.1/nu; % Reynolds number

%%
%--- set up ROM online parameters
snapIndex = 1; % initial condition, we choose first snapshot

T =10; % set integration time, [0,T]

dt=2.000e-03; % time step size % same as FOM

numTimeSteps = round(T/dt)-1; % number of time steps

Uave = PodU2d(:,1); % first 'mode' means the centering trajectory

% projection of the DNS snapshots on ROM basis
% u_proj = Mr\C*M_h*Y
dns_proj_data =MassROM(2:end,2:end)\( PodU2d(:,2:end)'*MassMatrix*(Snapshots(:,1:1500)-Uave));


r=8; % dimension of G-ROM

% now choose the first snapshot projected
uvproj0 = dns_proj_data(1:r,snapIndex);

[Gsolns,Mr] = nse_g_rom_online_fun(r,uvproj0,MassROM,StiffROM,TriLinRomD);


%%
%--- calculate kinetic energy

ke_g = zeros(1,numTimeSteps+1) ;
error_l2 = zeros(1,numTimeSteps+1) ;

tic
% calculate L2 error and kinetic energy

for i =1: numTimeSteps+1
    
    tmp_rom_vel =PodU2d(:,1) + PodU2d(:,2:r+1)*Gsolns(:,i);
    
    energy = 1/2 * sqrt(tmp_rom_vel' *MassMatrix * tmp_rom_vel );
    tmp_error = sqrt((tmp_rom_vel-Snapshots(:,i))' *MassMatrix * (tmp_rom_vel-Snapshots(:,i)) );
    
    ke_g(i) = energy;
    error_l2(i) = tmp_error;
end
toc

load('ke_dns_re1000.mat')

%%
t = 0:0.002:9.998;
bluecolor = [0 0.4470 0.7410];
redcolor = [0.8500 0.3250 0.0980];
%--- plot kinetic energy
figure(1); clf; jump = 50; frame_count = 1;
title('Kinetic Energy','FontSize',30,'FontName','times')
axis([0, 10, 0.6, 0.72]);
set(gcf, 'Position',  [0, 0, 2000, 550])
hold on
grid minor

for k = 1:jump:length(t)-jump
    plot(t(k:k+jump),ke_dns(k:k+jump),'Color','#A2142F','linewidth',4);
    hold on
    plot(t(k:k+jump),ke_g(k:k+jump),'Color','#0072BD','linewidth',4);
    ax = gca;
    ax.FontSize = 25;
    text(0.8, 1, ['t = ', num2str(t(k))]);
    M(frame_count) = getframe;
    hold on
    frame_count = frame_count + 1;
    
end
legend('FOM','G-ROM')


%--- plot l2 error

figure(2); clf; jump = 50; frame_count = 1;
title('L^2 Error','FontSize',30,'FontName','times')
axis([0, 10, 0, 0.7]);
set(gcf, 'Position',  [0, 650, 2000, 550])
hold on
grid minor
for k = 1:jump:length(t)-jump
    plot(t(k:k+jump),error_l2(k:k+jump),'Color','#7E2F8E','linewidth',4);
    ax = gca;
    ax.FontSize = 25;
    text(0.8, 1, ['t = ', num2str(t(k))]);
    M(frame_count) = getframe;
    hold on
    frame_count = frame_count + 1;
    
end






