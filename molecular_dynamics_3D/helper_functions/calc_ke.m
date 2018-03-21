%% calc_ke and 
%% kbT using equipartition
% why do we need the lattice constant in calculating ke?

function [kbT, kin_e] = calc_ke(velocities, lattice)
    v_squared = velocities.*velocities;
    kin_e = 0.5*sum(v_squared(:));
    [natoms, d] = size(velocities);
    kbT = (2/3)*(kin_e/natoms);%%divide by length to measure an average
end