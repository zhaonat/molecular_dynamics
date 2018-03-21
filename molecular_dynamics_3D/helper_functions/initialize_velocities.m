function [velocities,atoms_old] = initialize_velocities(atoms,latvec,kb_T,dt);
    %establish random velocities at temperature T

    [natoms,temp]=size(atoms);
    velocities=random('Uniform',-1,1,[natoms,3]);
    %velocities=0.1*ones(natoms, 3);

    %make center of mass velocity zero
    vel_cm=sum(velocities)/natoms;
    for k=1:natoms
        velocities(k,:)=velocities(k,:)-vel_cm;
    end

    kinetic_energy=0.5*sum(sum(velocities.*velocities));
    target_kinetic_energy=1.5*kb_T*natoms;
    scale_factor=sqrt(target_kinetic_energy/kinetic_energy);
    velocities=velocities*scale_factor;

    %set previous atom positions for use with original Verlet 
    %once we initialize the velocities, we technically know what the
    %positions of the atoms were -dt time previously
    %don't know why we need this?
    atoms_old=atoms-velocities*dt; %note that this is very different from the verlet integration
    
end





