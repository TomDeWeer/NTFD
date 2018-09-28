function fXi = calc_fXi(nF,cCoord,fNbCLoc,fNbC)
    for i = 1:nF
        leftCellIndex = fNbC(2*i-1);
        rightCellIndex = fNbC(2*i);
        leftCellCoord = fNbC(:,leftCellIndex);
        rightCellCoord = fNbC(:,rightCellIndex);
        Xi = rightCellCoord - leftCellCoord;
        % fnormal niet gegeven
    end
end



