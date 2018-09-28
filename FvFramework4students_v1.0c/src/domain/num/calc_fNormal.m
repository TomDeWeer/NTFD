function fNormal = calc_fNormal(fNbVLoc,fNbV,vCoord,fArea)
   % wijst van kleinste cell index naar grootste
    fNormal = [];
    for i = 1:nF
        vertexIndices = fNbV(2*i-1:2*i);
        firstLoc = vCoord(:,vertexIndices(1));
        secLoc = vCoord(:,vertexIndices(2));
        faceVector = secLoc-firstLoc;
        normalVector = [faceVector(2), -faceVector(1)];
        fNormal = [fNormal; normalVector];
    end
end