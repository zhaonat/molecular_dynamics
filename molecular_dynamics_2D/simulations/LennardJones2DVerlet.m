%% N particle, 2D verlet with lennard-jones potential
close all
clear;
NumParticles = 100;
Nx = 10;
Ny = 10;

%% Force Definition
% Fljx = @(rjk, xj, xk,r0) 1*(12*(r0/rjk)^13 - 6*(r0/rjk)^7)*(xj-xk)/rjk
% Fljy = @(rjk, yj, yk,r0) 1*(12*(r0/rjk)^13 - 6*(r0/rjk)^7)*(yj-yk)/rjk
% %% ====================
%% Plot of the potential


xVec = zeros(1,NumParticles); %% we'll update the positions without storing them
yVec = zeros(1, NumParticles);
counter = 1;
for i = 1:10
   
    for j = 1:10
       yVec(counter) = i;
       xVec(counter) = j; 
       counter = counter+1;  
    end
    
end
xinit = xVec; yinit = yVec;
% vxVec = randi([-10 10], NumParticles,1);
% vyVec = randi([-1 1], NumParticles,1);
vxVec = randn(NumParticles,1);
vyVec = randn(NumParticles,1);
vxinit = vxVec; vyinit = vyVec;
figure;
plot(xVec, yVec, '.', 'markersize', 30)
figure;
Nsteps = 5000;
tf = 50; t0 = 0;
r0 = 2; %radius where the lennard jones potential is 0
h = (tf-t0)/Nsteps;
eps = 1e-10;
for i = 1:Nsteps-1
    Ax1 = zeros(NumParticles,1);
    Ay1 = zeros(NumParticles,1);
    for j = 1:NumParticles 
       %% first half of the integration step (done for all particles)
       xj = xVec(j); yj = yVec(j);
       
       [ajx, ajy] = interparticleLJForce(xj, yj, xVec, yVec, NumParticles, eps, r0);
       Ax1(j) = ajx; Ay1(j) = ajy;
    end
    xold = xVec; yold = yVec;
    dr = zeros(NumParticles,2);
%     xVec = xVec+vxVec.'*h + 0.5*Ax1.'*h^2;
%     dr(:,1) = vxVec*h + 0.5*Ax1*h^2;
%     dr(:,2) = vyVec*h + 1/2*Ay1*h^2; 
    for j = 1:NumParticles
       xVec(j) = xVec(j) + vxVec(j)*h + 1/2*Ax1(j)*h^2;
       yVec(j) = yVec(j) +vyVec(j)*h + 1/2*Ay1(j)*h^2; 
       dr(j,1) = vxVec(j)*h + 1/2*Ax1(j)*h^2; 
       dr(j,2) = vyVec(j)*h + 1/2*Ay1(j)*h^2; 
    end
    Ax2 = zeros(NumParticles,1);
    Ay2 = zeros(NumParticles,1);
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
    
    if(mod(i, 2) == 0)
        plot(xVec,yVec, '.','markersize', 20)
        xlim([0 10])
        ylim([0,10])
        f = getframe;
    end
    i
end
figure;
plot(xVec, yVec, 'x')