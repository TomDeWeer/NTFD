%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%      NUMERICAL TECHNIQUES IN FLUID DYNAMICS       %%%%%%%%%%%%
%%%%%%%%%%%%             Koen Devesse, Tom De Weer             %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            SIMPLE algorithm                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% General test case
%
% A file used to quickly test the SIMPLE algorithm. This file should not be
% used to produce results: if figures or results for the report are needed,
% a proper test file for the specific case should be made.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

% Create a mesh
Lx = 1;
Ly = 1;
Nx = 20;
Ny = 20;
mu = 4;
Uxtop = 0;
p1 = 5;
p2 = 0;
dPdx = (p2-p1)/Lx;
rho = 10;
nu = mu/rho;
dx = Lx/Nx; 
dy = Ly/Ny;
seedI = LineSeed.lineSeedOneWayBias([0 0],[Lx 0],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 Ly],Ny,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');


% % actual solution
pfunc = @(x,y) p2 - dPdx*(Lx-x);
ufunc = @(x,y) (-dPdx/(2*mu))*y*(Ly-y)+Uxtop*y/Ly;
vfunc = @(x,y) 0;

% Set up initial fields
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[zeros(1,U.elcountzone);zeros(1,U.elcountzone)]);
casedef.U = U; % initial guess
P = Field(casedef.dom.allCells,0); % Pressure
P0 = [];
for i=1:casedef.dom.nC
    coord = casedef.dom.cCoord(:,i);
    x = coord(1);
    y = coord(2);
%     P0 = [P0, pfunc(x,y)];

%     P0 = [P0, pfunc(x,y) + sin(pi*x/dx )];
    P0 = [P0, sin(pi*x/dx )];
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
casedef.material.nu = nu;  % viscosity [kinematic]
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
casedef.BC{jBC}.data.velocity = [0, 0]; % uniform velocity inlet
casedef.BC{jBC}.pressureKind   = 'Dirichlet';
casedef.BC{jBC}.data.pressure = p1; % prescribed pressure
casedef.BC{jBC}.isNormalized = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.velocityKind   = 'Neumann';
casedef.BC{jBC}.data.velocity = [0, 0]; % developed velocity outlet, nonchanging
casedef.BC{jBC}.pressureKind   = 'Dirichlet';
casedef.BC{jBC}.data.pressure = p2; % prescribed pressure
casedef.BC{jBC}.isNormalized = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.velocityKind   = 'Dirichlet';
casedef.BC{jBC}.data.velocity = [0, 0]; % no slip condition
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.velocityKind   = 'Dirichlet';
casedef.BC{jBC}.data.velocity = [Uxtop, 0]; % no slip condition
casedef.BC{jBC}.pressureKind   = 'Neumann';
casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
casedef.BC{jBC}.isNormalized = 0;


% Set up iteration parameters
casedef.iteration.maxniter = 500;
casedef.iteration.resTol     = 1.e-5;
casedef.iteration.dt = 0.01;
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
for i=1:length(casedef.P.data)
    p = casedef.P.data(i);
    coord = casedef.dom.cCoord(:,i);
    x = coord(1);
    y = coord(2);
    if y == 0.5 + dx/2
        if x>0
            xi = [xi, x];
            pi = [pi, p];
        else
            xi = [x, xi];
            pi = [p, pi];
        end
    end
end
figure;
plot(xi, pi)

maxErr = 0;
avgErr = 0;
Ni=0;
Err = Field(casedef.dom.allCells, 0);
set(Err,zeros(1,U.elcountzone));
comp_err = zeros(1, size(Err.data,2));
ErrAtMiddle = [];
ValueAtMiddle = [];
ExactAtMiddle = [];
for i=1:result.U.dom.nC
    x = result.U.dom.cCoord(1,i);
    y = result.U.dom.cCoord(2,i);
    % Only keep the interior cells
    if x>0 && y>0 && x<Lx && y<Ly
        Uxapprox = result.U.data(1,i);
        Uxexact = ufunc(x,y);
        % Compute relative error
        err = abs((Uxexact-Uxapprox));
        comp_err(i) = err/abs(Uxexact);
        % Compute average error
        avgErr = avgErr + err;
        Ni = Ni+1;
        % Update maximum error
        if err>maxErr
            maxErr = err;
            maxErrx = x;
            maxErry = y;
        end
        if x == 0.5+dx/2
            ErrAtMiddle = [ErrAtMiddle, [y; err]];
            ValueAtMiddle = [ValueAtMiddle, [y; Uxapprox]];
            ExactAtMiddle = [ExactAtMiddle, [y; Uxexact]];
        end
    end
end
set(Err,comp_err);
figure; hold on; axis image; colormap(jet(50));
scale = 'lin'; lw = 0.2; colorbar(); title('Error')
fvmplotfield(Err,scale,lw);

avgErr = avgErr/Ni;
fprintf("Maximum error: %.10f \n", maxErr)
fprintf("Average error: %.10f \n", avgErr)

figure()
hold on
plot(ErrAtMiddle(1,:),ErrAtMiddle(2,:),'r');
plot(ValueAtMiddle(1,:),ValueAtMiddle(2,:),'b');
plot(ExactAtMiddle(1,:),ExactAtMiddle(2,:),'k');


