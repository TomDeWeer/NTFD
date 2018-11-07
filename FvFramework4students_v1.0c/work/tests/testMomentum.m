% First testcase: 1D flow from left to right, parabolic velocity inlet,
% converged velocity output, no pressure gradient???
clear all; close all;
% Create a mesh
seedI = LineSeed.lineSeedOneWayBias([0 0],[1 0],10,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 1],10,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');

% Set up initial fields
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[rand(1,U.elcountzone);rand(1,U.elcountzone)]);
% reset(U,[1;0.2]); 
% reset(U,[0; 0]);                          % Reset with all zeros
casedef.U = U; % initial guess
gradP = Field(casedef.dom.allCells,1); % Pressure
gradP0 = [];
for i=1:casedef.dom.nC
    % pos = casedef.dom.cCoord(i); % if its a function of position
    dPx = -1; % positive pressure gradient accelarating the flow
    dPy = 0;
    gradP0 = [gradP0, [dPx; dPy]];
end
set(gradP,gradP0)
casedef.gradP = gradP;
% Define material properties
casedef.material.nu = 1;  % viscosity [dynamic???]
casedef.material.rho = 1; % density [kg/m^3]


% Define boundary conditions
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
% casedef.BC{jBC}.kind   = 'Dirichlet';
% casedef.BC{jBC}.data.bcval = @(x,y) [y*(1-y), 0];
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = [0, 0];
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = [0 , 0];
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = [0, 0];
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = [0, 0];


% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.UTol     = 1e-6;
casedef.iteration.dt = 1.e-2;

% Call solver
result = momentumsolver(casedef);


% Plot result
% ux
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("Ux"); colorbar();
fvmplotfield(result.U,scale,lw, 1);
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("Uy"); colorbar();
fvmplotfield(result.U,scale,lw, 2);

figure()
quiver(casedef.dom.cCoord(1,:),...
    casedef.dom.cCoord(2,:),result.U.data(1,:),result.U.data(2,:));


