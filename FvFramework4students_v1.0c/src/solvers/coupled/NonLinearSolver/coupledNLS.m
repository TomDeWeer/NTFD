%==========================================================================
%
%
%
%==========================================================================
function sol = coupledNLS(casedef)
dom = casedef.dom;

% initialize the optimization vector x
% x contains velocities and pressures at every point in the grid
% x = [p ; u; v] where p contains pressures u contains x velocities and y
% contains y velocities

% uvec = casedef.U.data;
% u = uvec(1,:)';
% v = uvec(2,:)';
% p = casedef.P.data';
% 
% x0 = [p; u; v];
x0 = zeros(3*dom.nC,1);
outputFunc = @(x,optimVals,state) plotFlow(x,casedef, optimVals);
options = optimoptions('fsolve','Display','iter-detailed',...
    'OutputFcn',outputFunc,'SpecifyObjectiveGradient',true, ...
    'FiniteDifferenceStepSize', 1.e-5, 'FunctionTolerance',1.e-3, 'Algorithm','trust-region-dogleg'); % "CheckGradients",true,
handle = @(x) NavierStokes(casedef, x);
tic
sol = myOwn_fsolve(handle,x0,options);
fprintf("Time: %.3f \n",toc)
end
