% First testcase: square plate with expected peclet temperature decrease
% from left to right
clear all; close all;

path = [pwd '\Figuren\Convection'];
[~, ~, ~] = mkdir(path);
% Create a mesh
Nx = 150;
dx = 1/Nx;
L = 1;
V = 100;
D = 10; 
Pe = L*V/D;
fprintf("Peclet number: %.5f \n",Pe);
seedI = LineSeed.lineSeedOneWayBias([0 0],[L 0],Nx,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 1],1,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');

% Set up initial fields
T = Field(casedef.dom.allCells,0);     % Temperature [K] (scalar); empty field
% reset(T,0);                          % Reset with all zeros
randomdata = rand(T.elsize,T.elcountzone)-0.5;
set(T,randomdata);                     % Set with random numbers

U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[V*ones(1,U.elcountzone);zeros(1,U.elcountzone)]);
% reset(U,[1;0.2]); 
% reset(U,[0; 0]);                          % Reset with all zeros
casedef.U0 = U;

% Define material properties
casedef.material.k = D;  % Thermal conductivity [W/(m K)]

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
figure; hold on; axis image; colormap(jet(50));
scale = 'lin'; lw = 0.0;
fvmplotfield(result.T,scale,lw);
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
colorbar('TickLabelInterpreter', 'latex');
saveas(gcf,fullfile(path,'Numerical_results.png'));

% % Verifying accuracy
% % Getting the temperatures at the horizontal line x=0 to x=1
line = zeros(10,1); % contains temperatures at horizontal line
linexloc = zeros(10,1); % contains  x-coord at horizontal line
for i=1:result.T.dom.nC
    x = result.T.dom.cCoord(1,i);
    y = result.T.dom.cCoord(2,i);
    if x>0 && x<1 
        ix = round((x+(dx*0.5))/dx);
        line(ix) = result.T.data(i);
        linexloc(ix) = x;
    end
end
% calculate real velocity profile
% figure()
% plot(line, linexloc,'b')
hold on
realTempF = @(x) 1-((exp(Pe*x)-1)/(exp(Pe)-1));
realLine = zeros(size(line));
for i=1:length(line)
    realLine(i) = realTempF(linexloc(i));
end
% plot(realLine, linexloc,'r')
avgErr = mean(abs(realLine-line));
maxErr = max(abs(realLine-line));
fprintf("Maximum flow profile error: %.10f \n",maxErr)


AverageErrors = [];
MaximumErrors = [];
Nxs = [5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, ];
for Nx=Nxs
    dx = 1/Nx;
    L = 1;
    V = 100;
    D = 10; 
    Pe = L*V/D;
    fprintf("Peclet number: %.5f \n",Pe);
    fprintf("N_x: %f \n",Nx);
    seedI = LineSeed.lineSeedOneWayBias([0 0],[L 0],Nx,1.00,'o');
    seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 1],4,1.00,'o');
    casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
    mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
    % Create domain from mesh
    casedef.dom = newdomain(mesh,'MyDomain');

    % Set up initial fields
    T = Field(casedef.dom.allCells,0);     % Temperature [K] (scalar); empty field
    % reset(T,0);                          % Reset with all zeros
    randomdata = rand(T.elsize,T.elcountzone)-0.5;
    set(T,randomdata);                     % Set with random numbers

    U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
    set(U,[V*ones(1,U.elcountzone);zeros(1,U.elcountzone)]);
    % reset(U,[1;0.2]); 
    % reset(U,[0; 0]);                          % Reset with all zeros
    casedef.U0 = U;

    % Define material properties
    casedef.material.k = D;  % Thermal conductivity [W/(m K)]


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

    % Verifying accuracy
    % Getting the temperatures at the horizontal line x=0 to x=1, at y=0.125
    line = zeros(10,1); % contains temperatures at horizontal line
    linexloc = zeros(10,1); % contains  x-coord at horizontal line
    for i=1:result.T.dom.nC
        x = result.T.dom.cCoord(1,i);
        y = result.T.dom.cCoord(2,i);
        if y < 0.126 && y>0.124 && x>0 && x<1 % select if its at y=0.125
            ix = round((x+(dx*0.5))/dx);
            line(ix) = result.T.data(i);
            linexloc(ix) = x;
        end
    end
    % calculate real velocity profile
    realTempF = @(x) 1-((exp(Pe*x)-1)/(exp(Pe)-1));
    realLine = zeros(size(line));
    for i=1:length(line)
        realLine(i) = realTempF(linexloc(i));
    end
    avgErr = mean(abs(realLine-line));
    maxErr = max(abs(realLine-line));
    fprintf("Maximum flow profile error: %.10f \n",maxErr)
    fprintf("Average flow profile error: %.10f \n",avgErr)
    MaximumErrors = [MaximumErrors, maxErr];
    AverageErrors = [AverageErrors, avgErr];
end

figure()
loglog(Nxs,MaximumErrors, "kx")
hold on
loglog(Nxs, AverageErrors, "rx")
xlabel("$N_x$", 'Interpreter','latex')
ylabel("Error", 'Interpreter','latex')
legend("Maximum", "Average")
set(gca,'TickLabelInterpreter', 'latex');
saveas(gcf,fullfile(path,'Peclet_error_convergence.png'));


