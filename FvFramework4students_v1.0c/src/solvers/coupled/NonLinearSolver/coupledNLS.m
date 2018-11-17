%==========================================================================
%
%
%
%==========================================================================
function result = coupledNLS(casedef)
dom = casedef.dom;


% initialize the optimization vector x
% x contains velocities and pressures at every point in the grid
% x = [p ; u; v] where p contains pressures u contains x velocities and y
% contains y velocities
uvec = casedef.U.data;
u = uvec(1,:)';
v = uvec(2,:)';
p = casedef.P.data';

x0 = [p; u; v];

options = optimoptions('fsolve','SpecifyObjectiveGradient','on');
handle = @(x) NavierStokes(casedef, x);
sol = fsolve(handle,x0, options);



end
