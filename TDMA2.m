function T = TDMA2(T,aCoeff)

    [rows,cols] = size(T);
    
    %Defining source term Su
    Su = zeros(rows,cols); %No source
    
    %Extracting variables from struct variable
    aP = aCoeff.point;
    aE = aCoeff.east;
    aW = aCoeff.west;
    aN = aCoeff.north;
    aS = aCoeff.south;

    
    %Allocating memory for variables
    P = zeros(cols,1);
    Q = zeros(cols,1);
    
    for i = 2:rows-1
       
        %Treating special case for first finite element
        a = aP(i,2);
        b = aE(i,2);
        c = aW(i,2);
        d = aN(i,2)*T(i+1,2) + aS(i,2)*T(i-1,2) + Su(i,2);
        
        P(2) = b/a;
        Q(2) = (d+c*T(i,1))/a;
        
        for j = 3:cols-1 
            a = aP(i,j);
            b = aE(i,j);
            c = aW(i,j);
            d = aN(i,j)*T(i+1,j) + aS(i,j)*T(i-1,j) + Su(i,j);

            P(j) = b/(a - c*P(j-1));
            Q(j) = (d+c*Q(j-1))/(a - c*P(j-1));
        end
        
        %Solving system
        for j = cols-1:-1:2
            T(i,j) = Q(j) + P(j)*T(i,j+1);
        end
   
    end


end