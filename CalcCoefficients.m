function aCoeff = CalcCoefficients(T,x,y,u,v,rho,deltaX,deltaY,gamma,BC,...
    kFactor)

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
             
            %Calculating F and D
            Fw = (rho*u(i,j-1))*deltaY(i);
            Fe = (rho*u(i,j+1))*deltaY(i);
            Fs = (rho*v(i-1,j))*deltaX(j);  
            Fn = (rho*v(i+1,j))*deltaX(j);
            Dw  = kFactor*gamma*deltaY(i)/dXwest(j);
            De  = kFactor*gamma*deltaY(i)/dXeast(j);
            Ds  = kFactor*gamma*deltaX(j)/dYsouth(i);
            Dn  = kFactor*gamma*deltaX(j)/dYnorth(i);
            
            %Calculating coefficients
            west = max(Fw,(Dw+Fw/2));
            aCoeff.west(i,j) = max(0,west);
            east = max(-Fe,(De-Fe/2));
            aCoeff.east(i,j) = max(0,east);
            south = max(Fs,(Ds+Fs/2));
            aCoeff.south(i,j) =max(0,south);
            north = max(-Fn,(Dn-Fn/2));
            aCoeff.north(i,j) = max(0,north); 
          
        end
    end

    %Implementing boundary conditions in coefficients
    aCoeff.east(:,end-1) = aCoeff.east(:,end-1) * BC(2);
    aCoeff.west(:,2) = aCoeff.west(:,2) * BC(4);
    aCoeff.north(end-1,:) = aCoeff.north(end-1,:) * BC(3);
    aCoeff.south(2,:) = aCoeff.south(2,:) * BC(1);
    
    %Calculating aP
    aCoeff.point = aCoeff.east + aCoeff.west + aCoeff.south + aCoeff.north + Sp;

end