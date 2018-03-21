function [xnew, ynew, vxnew, vynew] = ReflectBC(rold, dr, vprev, Nx, Ny)
    x0 = rold(1);
    y0 = rold(2);
    dx = dr(1); dy = dr(2);
    vxnew = vprev(1); vynew = vprev(2);
    xnew = x0+dx;
    ynew = y0+dy;
    if(x0 +dx > Nx)
        diff = x0+dx-Nx;
        xnew = Nx-diff;
        vxnew = -vxnew;
    elseif ( x0+ dx < 0)
        xnew = 0+abs(x0+dx);
        vxnew = -vxnew;   
    end
    
    if(y0+dy > Ny)
        diff = y0+dy-Ny;
        ynew = Ny-diff;
        vynew = -vynew;
    elseif(y0+dy < 0)
        ynew = 0+abs(y0+dy);
        vynew = -vynew;   
    end
    

end