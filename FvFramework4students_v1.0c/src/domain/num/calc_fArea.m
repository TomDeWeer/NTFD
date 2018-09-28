function fArea = calc_fArea(nF,fNbVLoc,fNbV,vCoord)   
    fArea = [];
    for i = 1:nF
        vertexIndices = fNbV(2*i-1:2*i);
        firstLoc = vCoord(:,vertexIndices(1));
        secLoc = vCoord(:,vertexIndices(2));
        fAreai = norm(secLoc-firstLoc);
        fArea = [fArea; fAreai];
    end
    
end



