function [p, u, v] = getPUV(casedef, x)
nC = casedef.dom.nC;
p = x(1:nC);
u = x(nC+1:2*nC);
v = x(2*nC+1:end);

end