%% calculate correlation function essentially (should be roughly linear in time)
% inputs are atoms at time 0 and a list of saved atomic positions

function MSD = MSD(atoms_time_series, atoms0)
    [tsteps, numatoms, d] = size(atoms_time_series);
    MSD = [];
    for t = 1:tsteps
        disp = atoms_time_series(t,:,:)-atoms0;
        a = sum(sum(disp.*disp));
        MSD = [MSD, a];
    end

end