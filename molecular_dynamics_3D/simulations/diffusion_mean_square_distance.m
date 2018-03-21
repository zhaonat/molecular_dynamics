%% Problem 3 MD: STabilization of the Lattice TimeScale
%% Parameters L=M=N=3; T = 0.2 and T = 4; ONLY 20 TSTEPS!!!
addpath(genpath('D:\Nathan\Documents\StanfordYearOne\MatSci331\Homework3'))

%use LJ units
%atom positions are LJ units (not scaled)
close all
clear
dt = 0.01;
kb_T= 0.2;
nsteps=1000;

%set size of computational cell
L=2;
M=2;
N=2;
req = 2^(1/6);
lattice=sqrt(2)*2^(1/6);
%lattice=lattice*0.95;  %can scale lattice constant
rcut=3*req;
latvec=[L*lattice 0 0; 0 M*lattice 0; 0 0 N*lattice];

atoms = setup_cell(L,M,N, latvec);
[natoms, ~] = size(atoms);
%% PERFORM SIMULATION
[atomsf, instantaneous_kb_T,total_energy,pot_e, kin_e, saved_velocities,SD,...
    saved_trajectories] = ...
    runMD_MSD(kb_T, nsteps, L,M,N, rcut,dt);

%% visualize trajectories in 3D space
figure()
for i = 1:natoms
    scatter3(saved_trajectories(:,i,1), saved_trajectories(:,i,2), ...
        saved_trajectories(:,i,3), 'filled')
    hold on;
end
xlabel('x axis'); ylabel('y axis'), zlabel('z axis')
title(strcat('MD Trajectories at k_BT = ', num2str(kb_T), ...
    ' tsteps=',num2str(nsteps)))

%% Analysis
%SD = squared distance data;
[tsteps, natoms, dims] = size(SD);
MSD = [];

for i = 1:tsteps %% msd is a function of time
    value = 0;
    for d =1:dims %do the sum (xi-xi0)^2 for i = x,y,z
        %we sum over all the atoms and average
        value = value + sum(SD(i,:,d))/natoms;
    end
    MSD= [MSD, value];
    
end

figure()
plot(MSD, 'linewidth', 2)
grid()
title(strcat('Mean Square Displacement at k_bT=',num2str(kb_T)))
xlabel('time steps')
ylabel('MSD')
grid()

%% estimate slope using linear least squares
times = dt*(1:tsteps);
b1 = times.'\MSD.' %this assumes that the regression line goes through 0

figure()
%% matlab autocorr: time series for saved velocities
%integral of the autocorrelation gives another measure of D
for i = 1:natoms
    for j = 1:3
        mac = autocorr(saved_velocities(:,i,j),10);
        plot(mac)
        hold on;
    end
end
grid()


