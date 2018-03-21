function [atoms, instantaneous_kb_T,total_energy,pot_e, kin_e, saved_velocities] =...
    runMD(kb_T, nsteps, L,M,N, rcut,dt)
    %set size of computational cell

    %potential minimum is 2^(1/6)
    %set lattice constant (cubic primitive cell)
    lattice=sqrt(2)*2^(1/6);
    %lattice=lattice*0.95;  %can scale lattice constant

    %set lattice vectors
    latvec=[L*lattice 0 0; 0 M*lattice 0; 0 0 N*lattice];

    %set up computational cell for perfect xtal
    atoms=setup_cell(L,M,N,latvec);
    [natoms,temp]=size(atoms);

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
        [atoms,velocities, atoms_t]=integrate(atoms,rcut,velocities,forces,latvec,dt);
        
        %% MEASURE NEW KE AND PE
        if (rcut*2<latvec(1,1) && rcut*2<latvec(2,2) && rcut*2<latvec(3,3))
          [pot_e(time),forces]=calc_energy_faster(atoms,latvec,rcut,1);
        else
            [pot_e(time),forces]=calc_energy(atoms,latvec,rcut,1);
        end

        [instantaneous_kb_T(time),kin_e(time)]=calc_ke(velocities,lattice);
        total_energy(time)=kin_e(time)+pot_e(time);
        saved_velocities(time,:,:)=velocities;
        atoms_old=atoms;
        atoms=atoms;
        
        %% Plotting (optional)
%         s = 400;
%         scatter3(atoms(:,1), atoms(:,2), atoms(:,3),s,'filled')
%         xlim([0, L*lattice]); ylim([0,M*lattice]); zlim([0,N*lattice]);
%         getframe;
    end


end