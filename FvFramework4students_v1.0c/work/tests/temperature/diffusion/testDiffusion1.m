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
result = temperaturesolver(casedef);


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