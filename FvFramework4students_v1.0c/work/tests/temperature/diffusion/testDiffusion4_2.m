% Fourth testcase: influence of neumann boundary, sinusiodal flux
% at the top, isolated sides and grounded temperature at bottom
clear; close all; clc;

path = [pwd '/Figuren/Diffusion'];
[~, ~, ~] = mkdir(path);

H = 2;
L = 1;
number = 20;
K = 0.01;
Nx = number*L;
dx = L/Nx;
Ny = number*H;
dy = H/Ny;
% Create a mesh
seedI = LineSeed.lineSeedOneWayBias([0 0],[L 0],Nx,1,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 H],Ny,1,'o');
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
casedef.BC{jBC}.data.bcval = @(x,y) sin(pi*x/L); 
% de flux die door deze wand 
% naar binnen stroomt is 1W/m


% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.TTol     = 1e-6;


% Call solver
result = temperaturesolver(casedef);


% Plot result
figure; hold on; axis image; colormap(jet(50));
scale = 'lin'; lw = 0.2;
fvmplotfield(result.T,scale,lw);
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
colorbar('TickLabelInterpreter', 'latex');
saveas(gcf,fullfile(path,'Numerical_results.png'));

% Verifying accuracy
% analytic solution: eigenfunction
%An = @(n) 8*n/(K*(1-n^2)*(2*pi*n+1));
An = @(n) 4/(pi*K*(1-(n^2)));
A0 = 2*L/(pi*K);
% An = @(n) 2*sin(n*pi)/(K*n*pi);
% A0 = L/K;
% An = @(n) -4*(L^2)/(pi^2 * n^2 * K);
% A0 = L^3 / (6*K);
Tn = @(x,y,n) cos(n*pi*x/L).*sinh(n*pi*y/L)*An(n)*(L/(n*pi*cosh(n*pi*H/L)));
[X,Y] = meshgrid(dx/2:dx:L-dx/2,dy/2:dy:H-dy/2);
T = zeros(size(X));
T = T+A0*Y;
for n = 2:2:110
    %disp(An(n))
    Ti = Tn(X, Y, n);
    %disp(norm(Ti))
    if sum(isnan(Ti(:)))>0
       disp(n) 
       break
    end
    T = T+ Ti;
end
figure()
surf(X, Y, T)
shading interp
saveas(gcf,fullfile(path,'3D_surface.png'));

figure()
contourf(X, Y, T,20)
axis equal; colormap(jet(50));
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
colorbar('TickLabelInterpreter', 'latex');
saveas(gcf,fullfile(path,'Analytical_results.png'));

Ni=0;
maxErr = 0;
avgErr = 0;
TemperatureError = Field(casedef.dom.allCells, 0);
set(TemperatureError,zeros(1,TemperatureError.elcountzone));
comp_err = zeros(1, size(TemperatureError.data,2));
for i=1:result.T.dom.nC
    x = result.T.dom.cCoord(1,i);
    y = result.T.dom.cCoord(2,i);
    % Only keep the interior cells
    if x>0 && y>0 && x<L && y<H
        Tapprox = result.T.data(i);
        % Use knowledge from analyse 3 to compute exact solution
        Texact = A0*y;
        for n = 2:2:100
            Texact = Texact+ Tn(x, y, n);
        end
        % Compute relative error
        err = abs(Texact-Tapprox);
        comp_err(i) = err;

        % Compute average error
        avgErr = avgErr + err;
        Ni = Ni+1;
        % Update maximum error
        if err>maxErr
            maxErr = err;
            maxErrx = x;
            maxErry = y;
        end
    end
end
set(TemperatureError,comp_err);
figure; hold on; axis image; colormap(jet(50));
scale = 'lin'; lw = 0.2; colorbar(); title('Error')
fvmplotfield(TemperatureError,scale,lw);

avgErr = avgErr/Ni;
fprintf("Average error: %.10f \n",avgErr)
fprintf("Max error: %.10f \n",maxErr)
Ai = [];
ni = [];
for n = 2:2:110
    ni = [ni; n];
    Ai = [Ai; An(n)*max(max(Tn(X, Y, n)))];
end


figure()
semilogy(ni, abs(Ai),'b')
xlabel('n [-]','Interpreter','latex');
ylabel('T [K]','Interpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
saveas(gcf,fullfile(path,'Eigenfunction_size.png'));

%% Test convergentie

gridsizes = [2, 4, 8, 16, 32, 64, 128, 256];
error = [];
maxerrors = [];
for number=gridsizes
    disp(number)
    Nx = number*L;
    dx = L/Nx;
    Ny = number*H;
    dy = H/Ny;
    % Create a mesh
    seedI = LineSeed.lineSeedOneWayBias([0 0],[L 0],Nx,1,'o');
    seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 H],Ny,1,'o');
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
    casedef.BC{jBC}.data.bcval = @(x,y) sin(pi*x/L); 
    % de flux die door deze wand 
    % naar binnen stroomt is 1W/m


    % Set up iteration parameters
    casedef.iteration.maxniter = 1000;
    casedef.iteration.TTol     = 1e-6;


    % Call solver
    result = temperaturesolver(casedef);
    
    Ni=0;
    maxErr = 0;
    avgErr = 0;
    TemperatureError = Field(casedef.dom.allCells, 0);
    set(TemperatureError,zeros(1,TemperatureError.elcountzone));
    comp_err = zeros(1, size(TemperatureError.data,2));
    for i=1:result.T.dom.nC
        x = result.T.dom.cCoord(1,i);
        y = result.T.dom.cCoord(2,i);
        % Only keep the interior cells
        if x>0 && y>0 && x<L && y<H
            Tapprox = result.T.data(i);
            % Use knowledge from analyse 3 to compute exact solution
            Texact = A0*y;
            for n = 2:2:100
                Texact = Texact+ Tn(x, y, n);
            end
            % Compute relative error
            err = abs(Texact-Tapprox);
            comp_err(i) = err;
            % Compute average error
            avgErr = avgErr + err;
            Ni = Ni+1;
            % Update maximum error
            if err>maxErr
                maxErr = err;
                maxErrx = x;
                maxErry = y;
            end
        end
    end
    set(TemperatureError,comp_err);
    figure; hold on; axis image; colormap(jet(50));
    scale = 'lin'; lw = 0; colorbar();
    fvmplotfield(TemperatureError,scale,lw);
    
    figure; hold on; axis image; colormap(jet(50));
    scale = 'lin'; lw = 0;
    fvmplotfield(result.T,scale,lw);
    xlabel('x [m]','Interpreter','latex');
    ylabel('y [m]','Interpreter','latex');
    colorbar('TickLabelInterpreter', 'latex');
    saveas(gcf,fullfile(path,'Numerical_results.png'));

    avgErr = avgErr/Ni;
    maxerrors = [maxerrors, maxErr];
    error = [error, avgErr];
end

figure()
loglog(gridsizes, maxerrors, 'kx')
hold on
loglog(gridsizes, error, 'rx')
xlabel('N [-]','Interpreter','latex');
ylabel('Error [K]','Interpreter','latex');
h = legend('Maximum','Average');
set(h,'interpreter','Latex','FontSize',11)
set(gca,'TickLabelInterpreter', 'latex');
saveas(gcf,fullfile(path,'Error_convergence.png'));


















