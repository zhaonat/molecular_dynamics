%% Function to Calculate All InterParticle Forces

function [ajx,ajy] = interparticleLJForce(xj, yj, xVec, yVec, NumParticles, eps, r0)
     
    Fljx = @(rjk, xj, xk,r0, eps) eps*(12*(r0/rjk)^13 - 6*(r0/rjk)^7)*(xj-xk)/rjk;
    Fljy = @(rjk, yj, yk,r0, eps) eps*(12*(r0/rjk)^13 - 6*(r0/rjk)^7)*(yj-yk)/rjk;
    ajx = 0; ajy = 0;
    for k = 1:NumParticles
          xk = xVec(k); yk = yVec(k);
          %% xk == xj && yj ==yk is the case where we have the same particle, don't want that.
          if(xk == xj && yj == yk)
             continue; 
          end
          rjk = sqrt((xj-xk)^2 + (yj-yk)^2);
          ajx = ajx + Fljx(rjk, xj, xk,r0, eps);
          ajy = ajy + Fljy(rjk, yj, yk,r0, eps);

     end
     %Fxy = [ajx, ajy];
       
end