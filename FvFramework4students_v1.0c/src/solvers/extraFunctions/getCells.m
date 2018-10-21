function [ firstCell, secondCell ] = getCells( dom, indexFace )
%GETCELLS Returns the neighbouring cell indeces of a given face.
%
% Inputs:
%    dom    - Struct containing mesh information
%    indexFace  - The index of the given face
% 
% Outputs:
%    firstCell  - The index of the neighbouring cell with the lowest index.
%                 This will always be a physical cell.
%    secondCell - The index of the neighbouring cell with the highest
%                 index. If the given face is a boundary face, this will
%                 always be the ghost cell.
%

firstCell = dom.fNbC(2*indexFace-1);
secondCell = dom.fNbC(2*indexFace);

end

