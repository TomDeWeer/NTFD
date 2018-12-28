clear all; close all;
maxResiduals = [];
Ni = [20];
for N=Ni
    % Create a mesh
    Lx = 1;
    Ly = 1;
    Nx = N;
    Ny = N;
    mu = 4;
    dPdx = -10;
    Uxtop = 2;
    p0 = 5;
    rho = 10;
    nu = mu/rho;
    dx = Lx/Nx;
    dy = Ly/Ny;
    seedI = LineSeed.lineSeedOneWayBias([0 0],[Lx 0],Nx,1.00,'o');
    seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 Ly],Ny,1.00,'o');
    casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
    casedef.name = char("Channel__grid_"+Nx+"x"+Ny);

    mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
    % Create domain from mesh
    casedef.dom = newdomain(mesh,'MyDomain');

    % Define material properties
    casedef.material.nu = nu;  % viscosity [dynamic???]
    casedef.material.rho = rho; % density [kg/m^3]

    % Actual solution
    pfunc = @(x,y) p0 + dPdx*x;
    ufunc = @(x,y) (-dPdx/(2*mu))*y*(Ly-y)+Uxtop*y/Ly;
    vfunc = @(x,y) 0;

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
    casedef.BC{jBC}.data.velocity = [0,0]; % uniform velocity inlet
    casedef.BC{jBC}.pressureKind   = 'Dirichlet';
    casedef.BC{jBC}.data.pressure = pfunc; % prescribed pressure
    casedef.BC{jBC}.isNormalized = 0;
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'OOSTRAND';
    casedef.BC{jBC}.velocityKind   = 'Neumann';
    casedef.BC{jBC}.data.velocity = [0, 0]; % developed velocity outlet, nonchanging
    casedef.BC{jBC}.pressureKind   = 'Dirichlet';
    casedef.BC{jBC}.data.pressure = pfunc; % prescribed pressure
    casedef.BC{jBC}.isNormalized = 0;
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'ZUIDRAND';
    casedef.BC{jBC}.velocityKind   = 'Dirichlet';
    casedef.BC{jBC}.data.velocity = @(x,y) [ufunc(x,y), vfunc(x,y)]; % no slip condition
    casedef.BC{jBC}.pressureKind   = 'Neumann';
    casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
    casedef.BC{jBC}.isNormalized = 0;
    jBC = jBC+1;
    casedef.BC{jBC}.zoneID = 'NOORDRAND';
    casedef.BC{jBC}.velocityKind   = 'Dirichlet';
    casedef.BC{jBC}.data.velocity = @(x,y) [ufunc(x,y), vfunc(x,y)]; % no slip condition
    casedef.BC{jBC}.pressureKind   = 'Neumann';
    casedef.BC{jBC}.data.pressure = 0; % no normal pressure derivative
    casedef.BC{jBC}.isNormalized = 0;

    % Making real solution
    Ureal = Field(casedef.dom.allCells,1);
    data_Ureal = zeros(1, size(Ureal.data,2));
    Preal = Field(casedef.dom.allCells,0);
    data_Preal = zeros(1, size(Preal.data,2));

    for i=1:casedef.dom.nC
        x = casedef.dom.cCoord(1,i);
        y = casedef.dom.cCoord(2,i); 
        data_Preal(i) = pfunc(x,y);
        data_Ureal(2,i) = vfunc(x,y);
        data_Ureal(1,i) = ufunc(x,y);
    %     if (x>0 && y>0 && x<Lx && y<Ly) || (x<0) || (x>Lx)
    %         data_Ureal(1,i) = ufunc(x,y);
    %     elseif (y<0)
    %         data_Ureal(1,i) = -ufunc(x,-y);
    %     elseif (y>Ly)
    %         data_Ureal(1,i) = -ufunc(x,Ly-(y-Ly));
    %     end
    end
    set(Ureal,data_Ureal)
    set(Preal,data_Preal)
    p = data_Preal';
    u = data_Ureal(1,:)';
    v = data_Ureal(2,:)';
    res = NavierStokes(casedef,[p; u; v]);
    fprintf("Residual norm of analytical solution: %.10f \n",(res'*res))
    maxResiduals = [maxResiduals, max(abs(res))];
    figure; hold on; axis off; axis equal; colormap(jet(50));
    scale = 'lin'; lw = 0.2;
    set(gca,'TickLabelInterpreter', 'latex');
    xlabel('x [m]','Interpreter','latex');
    ylabel('y [m]','Interpreter','latex');
    colorbar('TickLabelInterpreter', 'latex');
    fvmplotfield(Preal,scale,lw);
end
% figure();
% loglog(Ni,maxResiduals, "kx")
% ylabel("max(residual)",'Interpreter','latex')
% set(gca,'TickLabelInterpreter', 'latex');