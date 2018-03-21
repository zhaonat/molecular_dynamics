%% N particle, 2D verlet with lennard-jones pot
%% Force Definition
% Fljx = @(rjk, xj, xk,r0) 1*(12*(r0/rjk)^13 - 6*(r0/rjk)^7)*(xj-xk)/rjk
% Fljy = @(rjk, yj, yk,r0) 1*(12*(r0/rjk)^13 - 6*(r0/rjk)^7)*(yj-yk)/rjk

function [atoms_tdt,velocities_tdt, atoms_old] = ...
    integrate(atoms, rcut, velocities, forces, latvec, dt);
    
    %% some simple initializations
    force_flag = 1;
    [NumParticles, d] = size(atoms);
    xVec = atoms(:,1); yVec = atoms(:,2); zVec = atoms(:,3);
    
    %% step 1: update particle positions
    atoms_old = atoms;
    atoms_tdt = atoms + velocities*dt + 0.5*forces*dt^2;
    
    %% step 2: get forces at t+dt using atoms_tdt;
    [E, forces_tdt] = calc_energy_faster(atoms_tdt,latvec, rcut, force_flag);

    %% step 2: update velocities to v_t+dt
    velocities_tdt = velocities + 0.5*(forces + forces_tdt)*dt;
 


    
end
