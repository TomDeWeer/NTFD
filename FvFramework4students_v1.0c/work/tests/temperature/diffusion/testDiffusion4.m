% Fourth testcase: influence of neumann boundary, sinusiodal flux
% at the top, isolated sides and grounded temperature at bottom
clear all; close all;
H = 0.5;
L = 1;

K = 0.01;
Nx = 50*L;
dx = L/Nx;
Ny = 50*H;
dy = H/Ny;
% Create a mesh
seedI = LineSeed.lineSeedOneWayBias([0 0],[L 0],Nx,0.95,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 H],Ny,0.9,'o');
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
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1;
fvmplotfield(result.T,scale,lw);

% Verifying accuracy
% analytic solution: eigenfunction
An = @(n) 8*n/(K*(1-n^2)*(2*pi*n+1));
A0 = 2*L/(pi*K);
Tn = @(x,y,n) cos(n*pi*x/L).*sinh(n*pi*y/L)*An(n)*(L/(n*pi*cosh(n*pi*H/L)));
[X,Y] = meshgrid(0:dx:L-dx,0:dy:H-dy);
X = X+dx/2;
Y = Y+dy/2;
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

figure()
contourf(X, Y, T,20)

Ni=0;
maxErr = 0;
avgErr = 0;
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
        err = abs(Texact-Tapprox)/Texact;
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
semilogy(ni, abs(Ai))

