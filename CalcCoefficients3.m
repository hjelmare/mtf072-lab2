function aCoeff = CalcCoefficients3(T,x,y,u,v,rho,deltaX,deltaY,gamma,BC,...
    kFactor,inletIndex,outletIndex)

    [rows,cols] = size(T);
    
    %Defining sourse term Sp
    Sp = zeros(rows,cols);  %No source
    
    %Allocating memory for variables
    aCoeff.east = zeros(rows,cols);
    aCoeff.west = zeros(rows,cols);
    aCoeff.north = zeros(rows,cols);
    aCoeff.south = zeros(rows,cols);
    deltaF = zeros(rows,cols);
    
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
            Fw = (rho*u(i,j-1))*deltaX(j)*deltaY(i)/dXwest(j);
            Fe = (rho*u(i,j+1))*deltaX(j)*deltaY(i)/dXeast(j);
            Fs = (rho*v(i-1,j))*deltaX(j)*deltaY(i)/dYsouth(i); 
            Fn = (rho*v(i+1,j))*deltaX(j)*deltaY(i)/dYnorth(i);
            deltaF(i,j) = Fe - Fw + Fn - Fs;
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

    %Implemenmting boundary conditions in coefficients
    aCoeff.east(:,end-1) = aCoeff.east(:,end-1) * BC(2);
    aCoeff.west(:,2) = aCoeff.west(:,2) * BC(4);
    aCoeff.north(end-1,:) = aCoeff.north(end-1,:) * BC(3);
    aCoeff.south(2,:) = aCoeff.south(2,:) * BC(1);
    aCoeff.point = aCoeff.east + aCoeff.west + aCoeff.south + ...
        aCoeff.north + Sp;
    
    %Implementing BC on inlet/outlet on west side
    j = 2;
    for i = inletIndex(1):inletIndex(end)
             
        %Calculating F and D
        Fw = (rho*u(i,j-1))*deltaX(j)*deltaY(i)/dXwest(j);
        Dw  = kFactor*gamma*deltaY(i)/dXwest(j);
        
        %Calculating coefficients
        aCoeff.west(i,j) = 2*Dw+Fw;  
        aCoeff.point(i,j) = aCoeff.east(i,j) + aCoeff.west(i,j) + ...
            aCoeff.south(i,j) + aCoeff.north(i,j) + Sp(i,j);
            
    end
   
    for i = outletIndex(1):outletIndex(end)
             
        %Calculating F and D
        Dw  = kFactor*gamma*deltaY(i)/dXwest(j);
        
        %Calculating coefficients
        aCoeff.west(i,j) = 2*Dw;  
        aCoeff.point(i,j) = aCoeff.east(i,j) + aCoeff.west(i,j) + ...
            aCoeff.south(i,j) + aCoeff.north(i,j) + Sp(i,j);
            
    end
end