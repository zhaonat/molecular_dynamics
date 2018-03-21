function [P_w,freqs] = calculate_Pw(saved_velocities,dt)
    %saved_velociites has dimension (iter x numPartilces x 3)
    [ntsteps, natoms, d] = size(saved_velocities);
    Fs = 1/dt;
    %% Step 1: fft the velocities for all particles
    P_w = 0;
    for ri = 1:d
        for ni = 1:natoms
           v_w = fftshift(fft(saved_velocities(:,ni,ri))); 
           P_w = P_w + abs(v_w).^2;
        end
    end
    P_w = P_w/(3*natoms*ntsteps); %normalize by time steps so the y scale is 1
    %% we need to generate the correct frequency axis for the spectrum as well
    freqs = Fs*(-(ntsteps/2):(ntsteps/2)-1)/ntsteps;
end