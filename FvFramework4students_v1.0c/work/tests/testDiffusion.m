% First testcase: square plate with expected linear temperature decrease
% from left to right
clear all; close all;
% Create a mesh
seedI = LineSeed.lineSeedOneWayBias([0 0],[1 0],10,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 1],10,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');

% Set up initial fields
T = Field(casedef.dom.allCells,0);     % Temperature [K] (scalar); empty field
% reset(T,0);                          % Reset with all zeros
randomdata = rand(T.elsize,T.elcountzone)-0.5;
set(T,randomdata);                     % Set with random numbers

U = Field(casedef.dom.allFaces,1);     % Velocity [m/s] (vector);
set(U,[rand(1,U.elcountzone);rand(1,U.elcountzone)]);
% reset(U,[1;0.2]); 
% reset(U,[0; 0]);                          % Reset with all zeros
casedef.U0 = U;

% Define material properties
casedef.material.k = 16;  % Thermal conductivity [W/(m K)]


% Define boundary conditions
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 1;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;


% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.TTol     = 1e-6;


% Call solver
result = examplesolver(casedef);


% Plot result
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1;
fvmplotfield(result.T,scale,lw);

% Verifying accuracy
% Getting the temperatures at the horizontal line x=0 to x=1, at y=0.55
line = zeros(10,1);
linexloc = zeros(10,1);
% Loop over all cells
for i=1:result.T.dom.nC
    x = result.T.dom.cCoord(1,i);
    y = result.T.dom.cCoord(2,i);
    % Select the cells at the horizontal line
    if y < 0.56 && y>0.54 && x>0 && x<1
        ix = round((x+0.05)/0.1);
        % Store temperature
        line(ix) = result.T.data(i);
        % Store cell x coordinate
        linexloc(ix) = x;
    end
end
% Plot results
figure()
plot(line, linexloc)
realTempF = @(x) (1-x); % Analytic solution(linear)
realLine = zeros(size(line));
for i=1:length(line)
    realLine(i) = realTempF(linexloc(i));
end
plot(realLine, linexloc)
err = norm(realLine-line)

% Second testcase: square plate isolated left and top, with T=1 at the
% right and T=0 at the bottom
clear all; close all;
% Create a mesh
N = 50
seedI = LineSeed.lineSeedOneWayBias([0 0],[1 0],N,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 1],N,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');

% Set up initial fields
T = Field(casedef.dom.allCells,0);     % Temperature [K] (scalar); empty field
% reset(T,0);                          % Reset with all zeros
randomdata = rand(T.elsize,T.elcountzone)-0.5;
set(T,randomdata);                     % Set with random numbers
% Define material properties
casedef.material.k = 16;  % Thermal conductivity [W/(m K)]
U = Field(casedef.dom.allFaces,1);     % Velocity [m/s] (vector);
set(U,[rand(1,U.elcountzone);rand(1,U.elcountzone)]);
casedef.U0 = U;
% Define boundary conditions
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 1;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.TTol     = 1e-6;
% Call solver
result = examplesolver(casedef);
% Plot result
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1;
fvmplotfield(result.T,scale,lw);

% Verifying accuracy
% analytic solution: eigenfunction
Tn = @(x,y,n) sin(pi*(0.5+n)*y).*cosh(pi*(0.5+n)*x)*2/(pi*(0.5+n)*cosh(pi*(0.5+n)));
[X,Y] = meshgrid(0:0.01:1,0:0.01:1);
T = zeros(size(X));
for n = 0:200
    T = T+ Tn(X, Y, n);
    if sum(isnan(T(:)))>0
       disp(n) 
    end
end
figure()
surf(X,Y,T)
figure()
contourf(X,Y,T)
% error calculation
% since the scheme is second order, the error should go down by factor 4 if
% dx and dy are simultaneously divided by 2
maxErr = 0;
avgErr = 0;
err0503 = 0;
% Loop over all cells
for i=1:result.T.dom.nC
    x = result.T.dom.cCoord(1,i);
    y = result.T.dom.cCoord(2,i);
    % Only keep the interior cells
    if x>0 && y>0 && x<1 && y<1
        Tapprox = result.T.data(i);
        % Use knowledge from analyse 3 to compute exact solution(first 200 terms)
        Texact = 0;
        for n = 0:200
            Texact = Texact+ Tn(x, y, n);
        end
        % Compute relative error
        err = abs(Texact-Tapprox)/Texact;
        % Compute average error
        avgErr = avgErr + err;
        % Update maximum error
        if err>maxErr
            maxErr = err;
            maxErrx = x;
            maxErry = y;
        end
        % Compute error at specific point(DOES NOT WORK YET)
        if x==0.5 && y==0.3
            err0503 = err
        end
    end
end
avgErr = avgErr/(N^2)
maxErr

% Third testcase: check the influence of changing cell sizes

% Fourth testcase: influence of neumann boundary, constant flux
% at the top, isolated sides and ground T at bottom
clear all; close all;
H = 2;
L = 1;
q = 1;
K = 17;
% Create a mesh
seedI = LineSeed.lineSeedOneWayBias([0 0],[L 0],10*L,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 H],100*H,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');

% Set up initial fields
T = Field(casedef.dom.allCells,0);     % Temperature [K] (scalar); empty field
% reset(T,0);                          % Reset with all zeros
randomdata = rand(T.elsize,T.elcountzone)-0.5;
set(T,randomdata);                     % Set with random numbers

U = Field(casedef.dom.allFaces,1);     % Velocity [m/s] (vector);
set(U,[zeros(1,U.elcountzone);zeros(1,U.elcountzone)]);
% reset(U,[1;0.2]); 
% reset(U,[0; 0]);                          % Reset with all zeros
casedef.U0 = U;

% Define material properties
casedef.material.k= K;  % Thermal conductivity [W/(m K)]


% Define boundary conditions
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = @(x,y) q; % de flux door deze wand is 1W/m


% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.TTol     = 1e-6;


% Call solver
result = examplesolver(casedef);


% Plot result
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1;
fvmplotfield(result.T,scale,lw);

% Verifying accuracy
% Getting the temperatures at the vertical line y=0 to y=1, at x=0.55
line = zeros(10,1);
lineyloc = zeros(10,1);
% Loop over all cells
for i=1:result.T.dom.nC
    x = result.T.dom.cCoord(1,i);
    y = result.T.dom.cCoord(2,i);
    % Select the cells at the horizontal line
    if x < 0.56 && x>0.54 && y>0 && y<=H
        ix = round((y+0.05)/0.1);
        % Store temperature
        line(ix) = result.T.data(i);
        % Store cell x coordinate
        lineyloc(ix) = y;
    end
end
% Plot results
figure()
hold on
plot(lineyloc, line)
realTempF = @(y) (q*y/K); % Analytic solution(linear), f(x)=1 so flux goes inside the domain
realLine = zeros(size(line));
for i=1:length(line)
    realLine(i) = realTempF(lineyloc(i));
end
plot(lineyloc, realLine)
err = norm(realLine-line)

% Fifth testcase: influence of neumann boundary, sinusiodal flux
% at the top, isolated sides and grounded temperature at bottom
clear all; close all;
H = 2;
L = 1;
q = 1;
K = 17;
Nx = 10*L;
dx = L/Nx;
Ny = 100*H;
dy = H/Ny;
% Create a mesh
seedI = LineSeed.lineSeedOneWayBias([0 0],[L 0],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 H],Ny,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');

% Set up initial fields
T = Field(casedef.dom.allCells,0);     % Temperature [K] (scalar); empty field
% reset(T,0);                          % Reset with all zeros
randomdata = rand(T.elsize,T.elcountzone)-0.5;
set(T,randomdata);                     % Set with random numbers

U = Field(casedef.dom.allFaces,1);     % Velocity [m/s] (vector);
set(U,[zeros(1,U.elcountzone);zeros(1,U.elcountzone)]);
% reset(U,[1;0.2]); 
% reset(U,[0; 0]);                          % Reset with all zeros
casedef.U0 = U;

% Define material properties
casedef.material.k= K;  % Thermal conductivity [W/(m K)]


% Define boundary conditions
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = @(x,y) sin(pi*x/L); % de flux door deze wand is 1W/m


% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.TTol     = 1e-6;


% Call solver
result = examplesolver(casedef);


% Plot result
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1;
fvmplotfield(result.T,scale,lw);

% Verifying accuracy
% analytic solution: eigenfunction
Tn = @(x,y,n) cos(n*pi*x/L).*sin(n*pi*y/L)*4*L*(cos(pi*n)+1)/((cos(n*pi*H/L))*(pi-pi*n^2)*(2*pi*n+sin(2*pi*n)));
[X,Y] = meshgrid(0:dx:L-dx,0:dy:H-dy);
X = X+dx/2;
Y = Y+dy/2;
T = zeros(size(X));
for n = 2:200
    T = T+ Tn(X, Y, n);
    if sum(isnan(T(:)))>0
       disp(n) 
    end
end
figure()
surf(X, Y, T)
disp('TODO: de juiste analytische oplossing vinden')