function T = GaussSeidel(T,aCoeff)

    [rows,cols] = size(T);
    
    %Defining source term Su
    Su = zeros(rows,cols); %No source
    
    %Extracting variables from struct variable
    aP = aCoeff.point;
    aE = aCoeff.east;
    aW = aCoeff.west;
    aN = aCoeff.north;
    aS = aCoeff.south;

    
    %Performing Gauss-Seidel update
    for j = 2:cols-1
        for i = 2:rows-1
            %Calculating parts of update equation
            Te = T(i,j+1);
            Tw = T(i,j-1);
            Tn = T(i+1,j);
            Ts = T(i-1,j);
            
            %Performing update on T(i,j)
            T(i,j) = (aE(i,j)*Te + aW(i,j)*Tw + aN(i,j)*Tn + aS(i,j)*Ts + Su(i,j))/aP(i,j);
                  
        end
    end

end