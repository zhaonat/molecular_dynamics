function [atoms, instantaneous_kb_T,total_energy,pot_e, kin_e, saved_velocities] =...
    runMD_thermostat(kb_T, nsteps, L,M,N, rcut,dt)
    %set size of computational cell
    kbTarget = kb_T;
    %potential minimum is 2^(1/6)
    %set lattice constant (cubic primitive cell)
    lattice=sqrt(2)*2^(1/6);
    %lattice=lattice*0.95;  %can scale lattice constant

    %set lattice vectors
    latvec=[L*lattice 0 0; 0 M*lattice 0; 0 0 N*lattice];

    %% set up computational cell for perfect xtal
    atoms=setup_cell(L,M,N,latvec);
    [natoms,temp]=size(atoms);
    eta = atoms;
    
    %% initialize velocities and measure initial values of everything
    [velocities,atoms_old]=initialize_velocities(atoms,latvec,kb_T,dt);
    %atoms_old is for time 0, so this is the missing velocity.
    [instantaneous_kb_T(1),kin_e(1)]=calc_ke(velocities,lattice);
    if (rcut*2<latvec(1,1) && rcut*2<latvec(2,2) && rcut*2<latvec(3,3))
    [pot_e(1),forces]=calc_energy_faster(atoms,latvec,rcut,1);
    else
    [pot_e(1),forces]=calc_energy(atoms,latvec,rcut,1);
    end
    total_energy(1)=kin_e(1)+pot_e(1);
    saved_velocities(1,:,:) = velocities;

    %% Enter Verlet Loop After Initialization Step
    Time_total = 0;
    for time=2:nsteps+1
        Time_total = Time_total+dt*time
        
        %% STEP EVERYTHING (atom positions and velocities) FORWARD
        [atoms,velocities, eta_t, atoms_t]=...
            integrate_thermostat(atoms,rcut,velocities,...
            latvec,dt, kbTarget, eta, lattice);
        eta = eta_t;
        %% MEASURE NEW KE AND PE
        if (rcut*2<latvec(1,1) && rcut*2<latvec(2,2) && rcut*2<latvec(3,3))
          [pot_e(time),forces]=calc_energy_faster(atoms,latvec,rcut,1);
        else
            [pot_e(time),forces]=calc_energy(atoms,latvec,rcut,1);
        end

        [instantaneous_kb_T(time),kin_e(time)]=calc_ke(velocities,lattice);
        total_energy(time)=kin_e(time)+pot_e(time);
        saved_velocities(time,:,:)=velocities;
                
        %% Plotting (optional)
        s = 400;
        if(mod(time, 10)==0)
            scatter3(atoms(:,1), atoms(:,2), atoms(:,3),s,'filled')
            xlim([0-1, L*lattice+1]); ylim([0-1,M*lattice+1]);
            zlim([0-1,N*lattice+1]);
            getframe;
        end
    end


end