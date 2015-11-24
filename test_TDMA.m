function [T] = test_TDMA(T,x,y,deltaX,deltaY,T1,c1,c2,kFactor)

    gamma = 1;
    
    [rows,cols] = size(T);  
    disp([rows,cols])
    
    %Predefining coefficients that handle boundary conditions
    eCoeff = ones(rows,cols);
    eCoeff(:,cols-1) = 2;
    wCoeff = ones(rows,cols);
    wCoeff(:,2) = 2;
    nCoeff = ones(rows,cols);
    nCoeff(rows-1,:) = 2;
    sCoeff = ones(rows,cols);
    sCoeff(2,:) = 2;
    
    % Create storage space
    a = zeros(cols,1);  % some of these can be shortened ...
    b = zeros(cols,1);  % and some are completely unecessary
    c = b;
    d = a;
    P = zeros(cols,1);
    Q = zeros(cols,1);
    
    for i = 2:rows-1 % for each row of the matrix storing the data points
        disp(['row ' num2str(i)])
        
        an = deltaY(i)/(y(i+1) - y(i)) * gamma;
        as = deltaY(i)/(y(i) - y(i-1)) * gamma;
        
        ae = deltaX(2)/(x(2+1) - x(2)) * gamma;
        aw = deltaX(2)/(x(2) - x(2-1)) * gamma;

        ap = ae + aw + an + as;
        
        a(2) = ap;
        b(2) = ae*eCoeff(i,2);
        c(2) = aw*wCoeff(i,2);
        d(2) = an*T(i-1,2) + as*T(i+1,2); % + source term??
                
        P(2) = b(2)/a(2);
        Q(2) = ( d(2)-c(2)*T(i,1) )/ a(2);
        for j = 3:cols    % for (almost) each column in the P-matrix
            ae = deltaX(j)/(x(j+1) - x(j)) * gamma;
            aw = deltaX(j)/(x(j) - x(j-1)) * gamma;
            
            ap = ae + aw + an + as;
            
            a(j) = ap;
            b(j) = ae*eCoeff(i,j);
            c(j) = aw*wCoeff(i,j);
            d(j) = an*T(i+1,j) + as*T(i-1,j); % + source term??
            
            P(j) = b(j) / (a(j) - c(j)*P(j-1));
            Q(j) = (d(j) + c(j)*Q(j-1))/(a(j)-c(j)*P(j-1));
        end
        
        %disp('abc')        
        %disp(diag(a(2:end),0) + diag(b(2:end-1),1) + diag(c(2:end-1),-1))
        %disp(' ')
        disp([a';b';c';d'])
        disp(' ')
        disp([P'; Q'])
        %disp(' ')
        
        T(i,end-1) = Q(end-1);
        for j = cols-2:-1:2
            %disp(['j ' num2str(j)]);
            T(i,j) = Q(j) + P(j)*T(i,j+1);
            %disp([Q(j), P(j), T(i,j+1)]);
        end
    end
    
    
end