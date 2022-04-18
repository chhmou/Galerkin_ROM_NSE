function [Gsolns, Mr] = nse_g_rom_online_fun(r,uvproj0,MassROM,StiffROM,TriLinRomD)
%{
This version does not have lift/drag/energy cal
This is Matlab function form 
%}

% \dot{a} = b + A*a + a^T*B*a; 

% Mr = MassROM(2:r+1,2:r+1);
% Sr = StiffROM(2:r+1,2:r+1);
% Tr = TriLinRomD(2:r+1,2:r+1,2:r+1);

%{
 VARIABLE NAME:
 Snapshots -- original UV field
 UVSnap -- UV perturbation field

 With given Mr, Sr, Gr and Tr
 Mr-- Mass matrix in ROM space, dim rxr
 Sr-- Stiffness matrix in ROM space, dim rxr
 Tr-- Trilinear term in ROM space, dim rxrxr
%}


global nu Re
global snapIndex
global T dt numTimeSteps

%format long e

Mr = MassROM(2:r+1,2:r+1);
Sr = StiffROM(2:r+1,2:r+1);
Tr = TriLinRomD(2:r+1,2:r+1,2:r+1);

Gsolns = zeros(r,numTimeSteps+1);
Gsolns(:,1) = uvproj0;
velPrevPrev = uvproj0;

ts  = 1;
%{
extract the (\Delta\phi_i,\Delta U) from Stiffness matrix
extract the (\phi_i,U\cdot U), (\phi_i,U\cdot \phi_j)
and (\phi_i,\phi_j\cdot U) from Trilinear term
%}
VecUs = StiffROM(2:r+1,1); % (\Delta\phi_i,\Delta U)


MatU1 = zeros(r,r);     % (\phi_i,U\cdot \phi_j)

for i =1:r
    temp1 = reshape(TriLinRomD(i+1,2:r+1,1),1,r);
    
    MatU1(i,:) = temp1;
    
end
VecUt0 = reshape(TriLinRomD(2:r+1,1,1),r,1); % (\phi_i,U\cdot\nabla U)

%%


VecUt = zeros(r,1); % (\phi_i,\phi_m\cdot\nabla U)
for i =1:r
    VecUt = VecUt+velPrevPrev(i)*reshape(TriLinRomD(2:r+1,1,i+1),r,1);
end

%%
%--------------------------------------------
% use the backward Euler to extropolate the u_r^{t1}
% velPrevPrev = velInit; %  u_r^{t0}

b = 1.0/dt * Mr * velPrevPrev-nu*VecUs-VecUt-VecUt0;


NLmat = 0*Mr;
for k=1:r
    NLmat = NLmat + velPrevPrev(k)*(Tr(:,:,k));
    %NLmat = NLmat + velPrevPrev(k)*(Tr(:,:,k)+ABtildeBnew(:,:,k)); % DD-VMS-ROM 
end
% linearized nonlinear term

A = 1.0/dt * Mr + nu*Sr + MatU1 + NLmat ;%+myAnew;
velSoln = A \ b;
velPrev =velSoln;
%velPrev= PodProjectionMatrix(:,snapIndex+1);% % test purpose

Gsolns(:,ts+1) = velSoln;



%% BDF2 with linear scheme

%
for ts=2:numTimeSteps
    %% Solve the ODE at each time step
    
    VecUt = zeros(r,1);
    
    for i =1:r
        VecUt = VecUt+(2*velPrev(i)-velPrevPrev(i))*reshape(TriLinRomD(2:r+1,1,i+1),r,1);
    end
    
    b = 2/dt * Mr * velPrev -0.5/dt*Mr*velPrevPrev -nu*VecUs-VecUt-VecUt0 ;
    % build matrix
    NLmat = 0*Mr;
    for k=1:r
        NLmat = NLmat + (2*velPrev(k)-velPrevPrev(k))*(Tr(:,:,k) );
        %NLmat = NLmat + (2*velPrev(k)-velPrevPrev(k))*(TriLinROM(:,:,k)+ABtildeBnew(:,:,k)); % DD-VMS-ROM 
    end
    A = 1.5/dt * Mr + nu*Sr + MatU1 + NLmat ;%+myAnew;
    % A = 1.5/dt * MassROM + nu*StiffROM + NLmat +myAnew; % DD-VMS-ROM 
    
    % solve the linear system
    velSoln = A \ b;
    
    
    Gsolns(:,ts+1)=velSoln;
    velPrevPrev=velPrev;
    velPrev = velSoln;
end

return

