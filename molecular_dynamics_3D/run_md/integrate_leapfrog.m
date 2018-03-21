%% N particle, 2D verlet with lennard-jones pot
%% Force Definition
% Fljx = @(rjk, xj, xk,r0) 1*(12*(r0/rjk)^13 - 6*(r0/rjk)^7)*(xj-xk)/rjk
% Fljy = @(rjk, yj, yk,r0) 1*(12*(r0/rjk)^13 - 6*(r0/rjk)^7)*(yj-yk)/rjk
%% Uses LeapFrog Verlet instead of velocity verlet, just as an alternate test
% implementation taken off of Wikipedia...
% note that the error in velocities is not quartic, but quadratic
function [atoms_tdt,velocities_tdt, atoms_old] = ...
    integrate_leapfrog(atoms, rcut, velocities, forces, latvec, dt);
    
    %% some simple initializations, information storage
    force_flag = 1;
    [NumParticles, d] = size(atoms);
    atoms_old = atoms;
    
    %% step 1: update positions by +1 step;
    atoms_tdt = atoms + velocities*dt;
    
    %% step 2: update forces
    [E, forces_dt] = calc_energy_faster(atoms_tdt,latvec, rcut, force_flag);

    %% step 2: update velocities to v_t+dt, by another half_step
    velocities_tdt = velocities + dt*forces_dt;
    
end
