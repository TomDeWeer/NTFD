clear all; close all;

% Create a mesh
Lx = 1;
Ly = 1;
Nx = 40;
Ny = 40;
mu = 4;
dPdx = -10;
Uxtop = 0;
p0 = 5;
rho = 10;
nu = mu/rho;
dx = Lx/Nx; 
dy = Ly/Ny;
theta = pi/4;
seedI = LineSeed.lineSeedOneWayBias([0 0],[sin(theta)*Lx cos(theta)*Lx],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[-cos(theta)*Ly sin(theta)*Ly],Ny,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');


% % actual solution
pfunc = @(x,y) p0 + dPdx*(-cos(theta)*x+sin(theta)*y);
ufunc = @(x,y) (-dPdx/(2*mu))*(cos(theta)*x+sin(theta)*y)*(Ly-(cos(theta)*x+sin(theta)*y))+Uxtop*(cos(theta)*x+sin(theta)*y)/Ly;
vfunc = @(x,y) 0;

% Set up initial fields
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[ones(1,U.elcountzone);ones(1,U.elcountzone)]);
casedef.U = U; % initial guess
P = Field(casedef.dom.allCells,0); % Pressure
P0 = [];
for i=1:casedef.dom.nC
    coord = casedef.dom.cCoord(:,i);
    x = coord(1);
    y = coord(2);
%     P0 = [P0, pfunc(x,y)];
    P0 = [P0, pfunc(x,y) + sin(pi*x/dx )];
%     P0 = [P0, 0];
end
set(P,P0)
casedef.P = P;

figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("P"); colorbar();
fvmplotfield(casedef.P,scale,lw);

xi = [];
pi = [];
for i=1:length(casedef.P.data)
    p = casedef.P.data(i);
    coord = casedef.dom.cCoord(:,i);
    x = coord(1);
    y = coord(2);
    if y<0.55 && y>0.45 && x>0 && x < Lx
        xi = [xi, x];
        pi = [pi, p];
    end
end
figure;
plot(xi, pi)


% Define material properties
casedef.material.nu = nu;  % viscosity [dynamic???]
casedef.material.rho = rho; % density [kg/m^3]


% Define boundary conditions
% There are boundary conditions for velocity and for pressure
% Wall: dirichlet velocity and neumann pressure in n
% Prescribed inlet: dirichlet velocity and dirichlet pressure
% Unprescribed outlet: neumann velocity in n and dirichlet pressure
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
% casedef.BC{jBC}.kind   = 'Dirichlet';
% casedef.BC{jBC}.data.bcval = @(x,y) [y*(1-y), 0];
casedef.BC{jBC}.velocityKind   = 'Neumann';
casedef.BC{jBC}.data.velocity = @(x,y) [0 , 0]; % uniform velocity inlet
casedef.BC{jBC}.pressureKind   = 'Dirichlet';
casedef.BC{jBC}.data.pressure = 5; % prescribed pressure
casedef.BC{jBC}.isNormalized = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.velocityKind   = 'Neumann';
casedef.BC{jBC}.data.velocity = [0, 0]; % developed velocity outlet, nonchanging
casedef.BC{jBC}.pressureKind   = 'Dirichlet';
casedef.BC{jBC}.data.pressure = 0; % prescribed pressure
casedef.BC{jBC}.isNormalized = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.velocityKind   = 'Dirichlet';
casedef.BC{jBC}.data.velocity = [0 , 0]; % no slip condition
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.velocityKind   = 'Dirichlet';
casedef.BC{jBC}.data.velocity = [0 , 0]; % no slip condition
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 0;


% Set up iteration parameters
casedef.iteration.maxniter = 150;
casedef.iteration.resTol     = 1.e-3;
casedef.iteration.dt = 20;
% relaxation factor
casedef.relaxation = 0.1;

% Call solver
tic
result = SIMPLEsolver(casedef);
fprintf("Time: %.3f \n",toc)

% Plot result
% ux
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("Ux"); colorbar();
fvmplotfield(result.U,scale,lw, 1);
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("Uy"); colorbar();
fvmplotfield(result.U,scale,lw, 2);
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("P"); colorbar();
fvmplotfield(result.P,scale,lw);

xi = [];
pi = [];
for i=1:length(result.P.data)
    p = result.P.data(i);
    coord = casedef.dom.cCoord(:,i);
    x = coord(1);
    y = coord(2);
    if y<0.55 && y>0.45 && x>0 && x < Lx
        xi = [xi, x];
        pi = [pi, p];
    end
end
figure;
plot(xi, pi)

