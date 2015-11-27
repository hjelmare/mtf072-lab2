function aCoeff = CalcCoefficients(T,x,y,deltaX,deltaY,gamma,BC,kFactor)

    [rows,cols] = size(T);
    
    %Defining sourse term Sp
    Sp = zeros(rows,cols);  %No source
    
    %Allocating memory for variables
    aCoeff.east = zeros(rows,cols);
    aCoeff.west = zeros(rows,cols);
    aCoeff.north = zeros(rows,cols);
    aCoeff.south = zeros(rows,cols);
    
    %Pre-calculating interpoint distances
    dXeast = diff(x(2:end));
    dXeast = [0 dXeast 0]';
    dXwest = diff(x(1:end-1));
    dXwest = [0 dXwest 0]';
    dYsouth = diff(y(1:end-1));
    dYsouth = [0 dYsouth 0]';
    dYnorth = diff(y(2:end));
    dYnorth = [0 dYnorth 0]';
    
    for j = 2:cols-1
        for i = 2:rows-1
            aCoeff.east(i,j) = kFactor * gamma * deltaY(i)/dXeast(j);
            aCoeff.west(i,j) = kFactor * gamma * deltaY(i)/dXwest(j);
            aCoeff.north(i,j) = kFactor * gamma * deltaX(j)/dYnorth(i);
            aCoeff.south(i,j) = kFactor * gamma * deltaX(j)/dYsouth(i);            
        end
    end

    %Implemenmting boundary conditions in coefficients
    aCoeff.east(:,2) = aCoeff.east(:,2) * BC(4);
    aCoeff.west(:,end-1) = aCoeff.west(:,end-1) * BC(2);
    aCoeff.north(end-1,:) = aCoeff.north(end-1,:) * BC(3);
    aCoeff.south(2,:) = aCoeff.south(2,:) * BC(1);
    
    %Calculating aP
    aCoeff.point = aCoeff.east + aCoeff.west + aCoeff.south + aCoeff.north + Sp;

end