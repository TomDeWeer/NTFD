function fXi = calc_fXi(nF,cCoord,fNbCLoc,fNbC)
    fXi = [];
    for i = 1:nF
        leftCellIndex = fNbC(2*i-1);
        rightCellIndex = fNbC(2*i);
        leftCellCoord = cCoord(:,leftCellIndex);
        rightCellCoord = cCoord(:,rightCellIndex);
        % moet liggen volgens normaal
        % dus van kleinste index naar grootste
        if leftCellIndex < rightCellIndex
            Xi = rightCellCoord - leftCellCoord;
        else
            Xi = leftCellCoord - rightCellCoord;
        end
        fXi = [fXi Xi];
    end
end



