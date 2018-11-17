
% Third testcase: influence of neumann boundary, constant flux
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
