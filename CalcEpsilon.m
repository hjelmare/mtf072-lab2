%function epsilon = CalcEpsilon(T,x,y,deltaX,deltaY,gamma,BC);

    [cols,rows] = size(T);
    
    %Source term
    bu = 0;
    bp = 0;
    
    %Precalculating aEast
    aEast = gamma .* (deltaY(2:end-1)./diff(x(1:end-1)));
    aEast(1) = aEast(1) * BC(4);
    
    %Precalculating aWest
    aWest = gamma .* (deltaY(2:end-1)./diff(x(2:end)));
    aWest(end) = aWest(end) * BC(2);
    
    %Precalculating aNorth
    aNorth = gamma .* (deltaX(2:end-1)./diff(y(2:end)));
    aNorth(end) = aNorth(end) * BC(3);
    
    %Precalculating aSouth
    aSouth = gamma .* (deltaX(2:end-1)./diff(y(1:end-1)));
    aSouth(1) = aSouth(1) * BC(1);
    
    for j = 2:cols-1
        ae = aEast(j);
        ew = aWest(j);
        for i = 2:rows-1
            Tp = T(i,j);
            Te = T(i,j+1);
            Tw = T(i,j-1);
            Tn = T(i+1,j);
            Ts = T(i-1,j);
            an = aNorth(i);
            as = aSouth(i);
            ap = ae + aw + an + as - bp;
            
            R = abs(ap*Tp - ae*Te - aw*Tw - an*Tn - as*Ts - bu);
            
        end
    end  
   





%end