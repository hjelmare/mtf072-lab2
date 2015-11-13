function [T] = initializeMesh(y,x,T1,T2)


    outMesh = zeros(length(y),length(x));
    outMesh(1,:) = T1;
    outMesh(:,end) = T2;
    outMesh(end,:) = 5;
    outMesh(:,1) = y*20;


    T = outMesh;
end


   