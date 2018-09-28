function fCentr = calc_fCentr(nF,fNbVLoc,fNbV,vCoord)
    fCentr = [];
    for i = 1:nF
        vertexIndices = fNbV(2*i-1:2*i);
        fCentri = zeros();
        for vertexIndex = vertexIndices
            vertexCoordinate = vCoord(:,vertexIndex);
            fCentri = fCentri + vertexCoordinate;
        end
        fCentri = fCentri / 2; % 2 is hardcoded
        fCentr = [fCentr; fCentri];
    end

end



