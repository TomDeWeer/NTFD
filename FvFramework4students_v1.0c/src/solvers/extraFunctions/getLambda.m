function [ lambda ] = getLambda( dom, indexFace )
%GETLAMBDA Determines lambda.
%   Used in at-face evaluation. The value at the face of variable Phi_face
%   can be linearly interpolated as:
%       Phi_face = lambda * Phi_cell1 + (1 - lambda) * Phi_cell2
%   This function determines and returns that lambda. Phi_cell1 is defined 
%   as the neighbouring cell with the lowest index.
% 
% Inputs:
%    dom    - Struct containing mesh information
%    indexFace  - The index of the given face
% 
% Outputs:
%    lambda - The value of lambda, as defined above.

[indexCell,~] = getCells(dom,indexFace);
lambda = 1 - norm(dom.cCoord(:,indexCell) - dom.fCoord(:,indexFace)) ...
    / dom.fXiMag(indexFace);

end

