function [ Au, bu, Av, bv ] = mMatrixMaker( casedef )
%MATRIXMAKER Returns the sparse matrix in vector form.
%   Detailed explanation goes here

dom = casedef.dom;

% Function to add the given terms to the appropriate places in the
% matrix.
function equationToUMatrix(diagIndex,offdiagIndex,diagTerm,offdiagTerm)
    Audiag(diagIndex) = Audiag(diagIndex) + diagTerm;    % Diagonal
    Auoffdiag(offdiagIndex) = offdiagTerm;   % Off-diagonal
end
function equationToVMatrix(diagIndex,offdiagIndex,diagTerm,offdiagTerm)
    Avdiag(diagIndex) = Avdiag(diagIndex) + diagTerm;    % Diagonal
    Avoffdiag(offdiagIndex) = offdiagTerm;   % Off-diagonal
end


nC = dom.nC;
nIf = dom.nIf;
nBf = dom.nBf;
Audiag = zeros(nC, 1);
Auoffdiag = zeros(2*(nIf+nBf), 1);
Avdiag = zeros(nC, 1);
Avoffdiag = zeros(2*(nIf+nBf), 1);
bu = zeros(nC,1);
bv = zeros(nC,1);
dt = casedef.iteration.dt;
rho = casedef.material.rho;
% add constant diagonal terms to A
for i=1:nC
   Ac = casedef.dom.cVol(i);% cell surface area
   Audiag(i) = Audiag(i) + Ac/dt;
   Avdiag(i) = Avdiag(i) + Ac/dt;
end
% add source terms
for i=1:nC
   Ac = casedef.dom.cVol(i);% cell surface area
   gradP = casedef.gradP.data(:,i);
   ui = casedef.U.data(1,i);
   vi = casedef.U.data(2,i);
   bu(i) = -(1/rho)*gradP(1)*Ac + Ac*(ui)/dt;
   bv(i) = -(1/rho)*gradP(2)*Ac + Ac*(vi)/dt;
end
% Compute coefficients for physical cell eqns and add them to eqn object
for i= 1:nIf+nBf
    % Getting terms of the equations
    [viscTerm,firstConvTerm,secondConvTerm] = MequationTerms(casedef,i);
    [firstCell,secondCell] = getCells(dom,i);
    % First cell
    firstDiag = viscTerm + firstConvTerm;
    firstOffdiag = - viscTerm - secondConvTerm;
    equationToUMatrix(firstCell,2*i-1,firstDiag,firstOffdiag);
    equationToVMatrix(firstCell,2*i-1,firstDiag,firstOffdiag);
    % Second cell
    % Checking whether it's a physical cell or a ghost cell
    if secondCell <= dom.nPc % Physical cell
        secondDiag = viscTerm + secondConvTerm;
        secondOffdiag = - viscTerm - firstConvTerm;
        equationToUMatrix(secondCell,2*i,secondDiag,secondOffdiag);
        equationToVMatrix(secondCell,2*i,secondDiag,secondOffdiag);
    else % If it's a ghost cell
        [ughostDiag, ughostOffdiag, ubValue, vghostDiag, vghostOffdiag, vbValue] ...
            = MghostTerms(casedef,i);
        equationToUMatrix(secondCell,2*i,ughostDiag,ughostOffdiag);
        equationToVMatrix(secondCell,2*i,vghostDiag,vghostOffdiag);
        bu(secondCell,:) = ubValue;
        bv(secondCell,:) = vbValue;
    end
end

Au = [Audiag; Auoffdiag];
Av = [Avdiag; Avoffdiag];

end