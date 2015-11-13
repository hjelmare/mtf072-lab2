function [T, pY, pX] = initializeMesh(edgesY, edgesX,T1,T2,T3,T4)


    edgesX = [2*edgesX(1)-edgesX(2) edgesX 2*edgesX(end)-edgesX(end-1)];
    edgesY = [2*edgesY(1)-edgesY(2) edgesY 2*edgesY(end)-edgesY(end-1)];

    ptPosX = edgesX(1:end-1) + 0.5*diff(edgesX);
    ptPosY = edgesY(1:end-1) + 0.5*diff(edgesY);

    outMesh = zeros(length(ptPosY),length(ptPosX));
    outMesh(1,:) = T1;
    outMesh(:,end) = T2;
    outMesh(end,:) = T3;
    outMesh(:,1) = T4;


    T = outMesh;
    pY = ptPosY;
    pX = ptPosX;

end


   