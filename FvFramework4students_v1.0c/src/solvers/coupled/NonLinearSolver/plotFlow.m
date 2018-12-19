function stop = plotFlow(x, casedef, optimVals)
[p, u, v] = getPUV(casedef,x);

Ux = Field(casedef.dom.allCells,0);     % Velocity [m/s] (x-comp);
set(Ux,u');
Uy = Field(casedef.dom.allCells,0);     % Velocity [m/s] (y-comp);
set(Uy,v');
P =  Field(casedef.dom.allCells,0); 
set(P,p');
Umag =  Field(casedef.dom.allCells,0); 
set(Umag,sqrt(u.^2+v.^2)');
outputpath = char("/users/start2015/r0585657/Documents/NTFD/NTFD/FvFramework4students_v1.0c/work/tests/momentum/coupled/NLS/Figures/LDC");

iteration = optimVals.iteration;
plot = 1;
if plot
    figure('visible','off'); hold on; axis off; axis equal; colormap(jet(50));
    scale = 'lin'; lw = 1; title("Ux"); colorbar();
    fvmplotfield(Ux,scale,lw, 1);
    saveas(gcf,fullfile(outputpath,char("Ux_"+iteration)),'png')
    figure('visible','off'); hold on; axis off; axis equal; colormap(jet(50));
    scale = 'lin'; lw = 1; title("Uy"); colorbar();
    fvmplotfield(Uy,scale,lw, 1);
    saveas(gcf,fullfile(outputpath,char("Uy_"+iteration)),"png")
    figure('visible','off'); hold on; axis off; axis equal; colormap(jet(50));
    scale = 'lin'; lw = 1; title("P"); colorbar();
    fvmplotfield(P,scale,lw, 1);
    saveas(gcf,fullfile(outputpath,char("P_"+iteration)),"png");
    close all;
    figure('visible','off'); hold on; axis off; axis equal; colormap(jet(50));
    scale = 'lin'; lw = 1; title("Umag"); colorbar();
    fvmplotfield(Umag,scale,lw, 1);
    saveas(gcf,fullfile(outputpath,char("Umag_"+iteration)),"png");
    close all;
    
    % plotting streamlines
%     dom = casedef.dom;
%     X = zeros();
%     Y = zeros();
%     UX = zeros();
%     UY = zeros();
%     for cellIndex = 1:dom.nC
%        pos = dom.cCoord(cellIndex);
%        
%     end
    
    
%     [x,y] = meshgrid(0:0.1:1,0:0.1:1);
%     u = x;
%     v = -y;
%     figure
%     quiver(x,y,u,v)
% 
%     startx = 0.1:0.1:1;
%     starty = ones(size(startx));
%     streamline(x,y,u,v,startx,starty)
end

stop=0; % choose stopping criterion here
end