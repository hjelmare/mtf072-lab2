function epsilon = CalcEpsilon(T,aCoeff,y)

    h = 0.03;   % Height of inlet/outlet

    [rows,cols] = size(T);
    
    %Defining source term Su
    Su = zeros(rows,cols); %No source
    
    %Extracting variables from struct variable
    aP = aCoeff.point;
    aE = aCoeff.east;
    aW = aCoeff.west;
    aN = aCoeff.north;
    aS = aCoeff.south;

    
    %Summing up residuals forall cells
    sumR = 0;
    for j = 2:cols-1
        for i = 2:rows-1
            %Extracting
            Te = T(i,j+1);
            Tw = T(i,j-1);
            Tn = T(i+1,j);
            Ts = T(i-1,j);
            
            %Calculating residual for cell i,j
            R = abs(aP(i,j)*T(i,j) - (aE(i,j)*Te + aW(i,j)*Tw + ...
                aN(i,j)*Tn + aS(i,j)*Ts + Su(i,j)));
                  
            sumR = sumR + R;
        end
    end

    flux = h * abs(mean(T(y<h,2)) - mean(T(y>2-h,2)));
    epsilon = sumR/flux;

end