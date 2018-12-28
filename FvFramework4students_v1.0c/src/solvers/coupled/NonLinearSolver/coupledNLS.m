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

% uvec = casedef.U.data;
% u = uvec(1,:)';
% v = uvec(2,:)';
% p = casedef.P.data';
% 
% x0 = [p; u; v];
x0 = zeros(3*dom.nC,1);
%outputFunc = @(x,optimVals,state) plotFlow(x,casedef, optimVals);
options = optimoptions('fsolve','Display','iter-detailed',...
    'SpecifyObjectiveGradient',true, ...
    'FiniteDifferenceStepSize', 1.e-5, 'FunctionTolerance',casedef.iteration.restol, 'Algorithm','trust-region-dogleg'); %"CheckGradients",true); % 'OutputFcn',outputFunc,
handle = @(x) NavierStokes(casedef, x);
tic
[sol, fval, exitflag, output] = myOwn_fsolve(handle,x0,options);
output.time = toc;
plotFlowFinal(sol,casedef,output)
[p, u, v] = getPUV(casedef,sol);
U = Field(casedef.dom.allCells,1);
set(U, [u' ; v']);
result.U = U;
P = Field(casedef.dom.allCells,0);
set(P, p')
result.P = P;
end
