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

% U = Field(casedef.dom.allFaces,1);     % Velocity [m/s] (vector);
% set(U,[rand(1,U.elcountzone);rand(1,U.elcountzone)]);
% reset(U,[1;0.2]); 
% reset(U,[0; 0]);                          % Reset with all zeros


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
for i=1:result.T.dom.nC
    x = result.T.dom.cCoord(1,i);
    y = result.T.dom.cCoord(2,i);
    if y < 0.56 && y>0.54 && x>0 && x<1
        ix = round((x+0.05)/0.1);
        line(ix) = result.T.data(i);
        linexloc(ix) = x;
    end
end
figure()
plot(line, linexloc)
realTempF = @(x) (1-x);
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
N = 200
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
% analytic solution
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
for i=1:result.T.dom.nC
    x = result.T.dom.cCoord(1,i);
    y = result.T.dom.cCoord(2,i);
    if x>0 && y>0 && x<1 && y<1
        Tapprox = result.T.data(i);
        Texact = 0;
        for n = 0:200
            Texact = Texact+ Tn(x, y, n);
        end
        err = abs(Texact-Tapprox)/Texact;
        avgErr = avgErr + err;
        if err>maxErr
            maxErr = err;
            maxErrx = x;
            maxErry = y;
        end
        if x==0.5 && y==0.3
            err0503 = err
        end
    end
end
avgErr = avgErr/(N^2)
maxErr

