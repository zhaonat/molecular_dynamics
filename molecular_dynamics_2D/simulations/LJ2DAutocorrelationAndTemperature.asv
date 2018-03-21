%% N particle, 2D verlet with lennard-jones potential
close all
clear;
NumParticles = 64;
Nx = 8;
Ny = 8;
Nsteps = 3000;
tf = 30; t0 = 0;
r0 = 2; %radius where the lennard jones potential is 0
h = (tf-t0)/Nsteps;
eps = 1e-10;

%% Force Definition
% Fljx = @(rjk, xj, xk,r0) 1*(12*(r0/rjk)^13 - 6*(r0/rjk)^7)*(xj-xk)/rjk
% Fljy = @(rjk, yj, yk,r0) 1*(12*(r0/rjk)^13 - 6*(r0/rjk)^7)*(yj-yk)/rjk
% %% ====================
%% Plot of the potential
AVG = [];

for T = 0.05:0.02:2.0

        xVec = zeros(1,NumParticles); %% we'll update the positions without storing them
        yVec = zeros(1, NumParticles);
        counter = 1;
        xTrack = zeros(NumParticles,Nsteps); yTrack = zeros(6,Nsteps);
        Autocorr = zeros(NumParticles,Nsteps);
        for i = 1:8

            for j = 1:8
               yVec(counter) = i;
               xVec(counter) = j; 

               xTrack(counter,1) = i;
               yTrack(counter,1) = j;

               counter = counter+1;
            end

        end
        xinit = xVec; yinit = yVec;
        % vxVec = randi([-10 10], NumParticles,1);
        % vyVec = randi([-1 1], NumParticles,1);
        vxVec = T*randn(NumParticles,1);
        vyVec = T*randn(NumParticles,1);
        vxinit = vxVec; vyinit = vyVec;

        %% calculate temperature using rootmean square
        vxMeanSq = mean(vxinit.^2);
        vyMeanSq = mean(vyinit.^2);
        vxrms = sqrt(vxMeanSq); vyrms = sqrt(vxMeanSq);
        vrms = sqrt(vxrms^2 + vyrms^2);
        T = vrms^2/3
        %%==========================================


        figure;

        %% full data structure
        x_t = zeros(NumParticles, Nsteps);

        %% Simulation Begins
        for i = 1:Nsteps-1
            Ax1 = zeros(NumParticles,1);
            Ay1 = zeros(NumParticles,2);
            for j = 1:NumParticles 
               %% first half of the integration step (done for all particles)
               xj = xVec(j); yj = yVec(j);

               [ajx, ajy] = interparticleLJForce(xj, yj, xVec, yVec, NumParticles, eps, r0);
               Ax1(j) = ajx; Ay1(j) = ajy;
            end
            xold = xVec; yold = yVec;
            dr = zeros(NumParticles,2);
            track = 1;
            for j = 1:NumParticles
               xVec(j) = xVec(j) + vxVec(j)*h + 1/2*Ax1(j)*h^2;
               yVec(j) = yVec(j) +vyVec(j)*h + 1/2*Ay1(j)*h^2; 
               dr(j,1) = vxVec(j)*h + 1/2*Ax1(j)*h^2; 
               dr(j,2) = vyVec(j)*h + 1/2*Ay1(j)*h^2; 

               %% store tracking variables and Extract Autocorrelation
               xTrack(j,i) = xVec(j);
               yTrack(j,i) = yVec(j);

               for k = 1:NumParticles
                  autoc =  AutoCorrBrownian(xTrack(k,:), yTrack(k,:), i);
                  Autocorr(k,i) = autoc;
               end

            end
            Ax2 = zeros(NumParticles,1);
            Ay2 = zeros(NumParticles,2);
            for(j= 1:NumParticles)
               %% second half of the step (done for all particles)
               xj = xVec(j); yj = yVec(j);

               [ajx2, ajy2] = interparticleLJForce(xj, yj, xVec, yVec, NumParticles, eps, r0);
               Ax2(j) = ajx2; Ay2(j) = ajy2;


            end
            for j = 1:NumParticles
               vxVec(j) = vxVec(j) + 1/2*(Ax1(j) + Ax2(j))*h;
               vyVec(j) = vyVec(j) + 1/2*(Ay1(j)+ Ay2(j))*h;

               %% ENFORCE REFLECTING BC
               [xnew, ynew, vxnew, vynew] = ReflectBC([xold(j), yold(j)], dr(j,:), [vxVec(j), vyVec(j)], Nx, Ny);
               xVec(j) = xnew;
               yVec(j) = ynew;
               vxVec(j) = vxnew;
               vyVec(j) = vynew;
            end    

            if(mod(i, 10) == 0)
                plot(xVec,yVec, '.','markersize', 25)
                xlim([0 8])
                ylim([0 8])
                xlabel('x coordinate')
                ylabel('y coordinate')
                title('MD Lennard-Jones Simulation')
                mov(i/10) = getframe;

            end
            x_t(:,i) = xVec;
            i
        end

        figure;
        plot(xVec, yVec, 'x')

        %% Autocorellation results
        k = figure;
        for i = 1:5
           plot(xTrack(i, :), yTrack(i,:), '.', 'markersize', 10)
           hold on;
        end
        xlabel('x coordinate')
        ylabel('y coordinate')
        title('Trajectories of Several Particles in MD Simulation')
        saveas(k, strcat('Trajectories for tf =', num2str(tf),' and T=', num2str(round(100*T))));

        Time = 0:h:tf-h;
        jf = figure;
        plot(Time, Autocorr.')
        xlabel('time')
        ylabel('Autocorrelation Value')
        %% Get Average
        avg = zeros(Nsteps,1);
        for i = 1:Nsteps
            val = 0;
            for j = 1:NumParticles
               val = val + Autocorr(j,i); 
            end
            val = val/NumParticles;
            avg(i) = val;
        end

        hold on;
        plot(Time, avg, '.-', 'markersize', 15);
        ylim([0 4]);
        title(strcat('Autocorrelation Plt for MD of LJ Potential at T=',num2str(round(100*T))))
        saveas(jf, strcat('Autocorrelation Plt at T=', num2str(round(100*T))))
        AVG = [AVG,avg];
        

end
figure;
plot(AVG, '.-')
title('Autocorrelation Plots at Different Temperatures')
xlabel('time')
ylabel('value')