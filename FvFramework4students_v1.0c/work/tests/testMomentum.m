% First testcase: 1D flow from left to right, parabolic velocity inlet,
% converged velocity output, no pressure gradient???
clear all; close all;

% Create a mesh
Lx = 1;
Ly = 1;
R = Ly/2;
Nx = 20;
Ny = 20;
nu = 1;
rho = 1;
dx = Lx/Nx;
dy = Ly/Ny;
seedI = LineSeed.lineSeedOneWayBias([0 0],[Lx 0],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 Ly],Ny,1.00,'o');
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
casedef.material.nu = nu;  % viscosity [dynamic???]
casedef.material.rho = rho; % density [kg/m^3]


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
casedef.iteration.UTol     = 1.e-10;
casedef.iteration.dt = 2;

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

% comparing to analytical solution
line = zeros(Ny,1);
lineyloc = zeros(Ny,1);
% Loop over all cells
for i=1:result.U.dom.nC
    x = result.U.dom.cCoord(1,i);
    y = result.U.dom.cCoord(2,i);
    % Select the cells at the vertical line
    if x < (0.5*Lx+dx) && x>(0.5*Lx) && y>0 && y<Ly
        iy = round((y+0.5*dy)/dy);
        % Store temperature
        line(iy) = result.U.data(1,i);
        % Store cell x coordinate
        lineyloc(iy) = y;
    end
end
% Plot results
figure()
hold on;
plot(lineyloc, line)
xlabel("y")
ylabel("u_x")
title("Fully developed flow profile")

realU = @(y) (R^2/(4*nu))*-dPx*(1-((R-y)/R)^2); % Analytic solution(linear)
realLine = zeros(size(line));
for i=1:length(line)
    realLine(i) = realU(lineyloc(i));
end
plot(lineyloc, realLine)
err = norm(realLine-line)



