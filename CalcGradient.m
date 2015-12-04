function [dX,dY] = CalcGradient(T,x,y)

    %Modifying input variables
    T = T(2:end-1,2:end-1);
    x = x(2:end-1);
    y = y(2:end-1);

    [rows,cols] = size(T);

    %Estimate gradient in every point with finite differences
    dX = zeros(rows,cols);
    dY = zeros(rows,cols);
    for j = 2:cols-1
        for i = 1:rows
            dX(i,j) = (T(i,j+1) - T(i,j-1)) / (x(j+1) - x(j-1));
        end
    end
    
    for j = 1:cols
        for i = 2:rows-1
            dY(i,j) = (T(i+1,j) - T(i-1,j)) / (y(i+1) - y(i-1));
        end
    end

end