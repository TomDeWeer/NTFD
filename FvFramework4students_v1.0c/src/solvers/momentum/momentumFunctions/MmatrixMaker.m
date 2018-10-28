function [ A, b ] = matrixMaker( casedef )
%MATRIXMAKER Returns the sparse matrix in vector form.
%   Detailed explanation goes here

dom = casedef.dom;

    % Function to add the given terms to the appropriate places in the
    % matrix.
    function equationToMatrix(diagIndex,offdiagIndex,diagTerm,offdiagTerm)
        Adiag(diagIndex) = Adiag(diagIndex) + diagTerm;    % Diagonal
        Aoffdiag(offdiagIndex) = offdiagTerm;   % Off-diagonal
    end

nC = dom.nC;
nIf = dom.nIf;
nBf = dom.nBf;
Adiag = zeros(nC, 1);
Aoffdiag = zeros(2*(nIf+nBf), 1);
b = zeros(nC,1);
% Compute coefficients for physical cell eqns and add them to eqn object
for i= 1:nIf+nBf
    % Getting terms of the equations
    [anb,firstConvTerm,secondConvTerm] = equationTerms(casedef,i);
    [firstCell,secondCell] = getCells(dom,i);
    % First cell
    firstDiag = anb + firstConvTerm;
    firstOffdiag = - anb - secondConvTerm;
    equationToMatrix(firstCell,2*i-1,firstDiag,firstOffdiag);
    % Second cell
    % Checking whether it's a physical cell or a ghost cell
    if secondCell <= dom.nPc % Physical cell
        secondDiag = anb + secondConvTerm;
        secondOffdiag = - anb - firstConvTerm;
    else % If it's a ghost cell
        [secondDiag,secondOffdiag,bValue] ...
            = ghostTerms(casedef,i);
        b(secondCell) = bValue;
    end
    % Filling in the offdiagonal elements
    equationToMatrix(secondCell,2*i,secondDiag,secondOffdiag);
end

A = [Adiag; Aoffdiag];

end

