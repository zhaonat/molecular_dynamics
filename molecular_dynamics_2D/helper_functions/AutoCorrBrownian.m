function autoc = AutoCorrBrownian(xVec, yVec, N)
    
    x0 = xVec(1);
    y0 = yVec(1);
    rsq = 0;
    for i = 2:N
      rsq = rsq +(xVec(i)-x0)^2 + (yVec(i) -y0)^2;
    end
    autoc = rsq/N;
    
end