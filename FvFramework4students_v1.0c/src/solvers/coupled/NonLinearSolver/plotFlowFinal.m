function stop = plotFlowFinal(x, casedef, output)
[p, u, v] = getPUV(casedef,x);

Ux = Field(casedef.dom.allCells,0);     % Velocity [m/s] (x-comp);
set(Ux,u');
Uy = Field(casedef.dom.allCells,0);     % Velocity [m/s] (y-comp);
set(Uy,v');
P =  Field(casedef.dom.allCells,0); 
set(P,p');
Umag =  Field(casedef.dom.allCells,0); 
set(Umag,sqrt(u.^2+v.^2)');

%outputpath = char("/users/start2015/r0585657/Documents/NTFD/NTFD/FvFramework4students_v1.0c/work/tests/momentum/coupled/NLS/Figures/LDC/"+casedef.name);
outputpath = char("C:\Users\Tom\Desktop\Universiteit\2emaster\NTFD\FvFramework4students_v1.0c\FvFramework4students_v1.0c\work\tests\momentum\coupled\NLS\Figures\Channel\"+casedef.name);
mkdir(outputpath)
save(char(outputpath+"/output.mat"),'output','x')
visible = "off";
plot = 1;
if plot
    figure('visible',visible); hold on; axis off; axis equal; colormap(jet(50));
    scale = 'lin'; lw = 0; title("Ux"); colorbar(); grid off; shading interp;
    fvmplotfield(Ux,scale,lw, 1);
%     saveas(gcf,fullfile(outputpath,'Ux.png'))
    figure('visible',visible); hold on; axis off; axis equal; colormap(jet(50));
    scale = 'lin'; lw = 0; title("Uy"); colorbar(); grid off; shading interp;
    fvmplotfield(Uy,scale,lw, 1);
%     saveas(gcf,fullfile(outputpath,'Uy.png'))
    figure('visible',visible); hold on; axis off; axis equal; colormap(jet(50));
    scale = 'lin'; lw = 0; title("P"); colorbar(); grid off; shading interp;
    fvmplotfield(P,scale,lw, 1);
%     saveas(gcf,fullfile(outputpath,'P.png'));
    %close all;
    figure('visible',visible); hold on; axis off; axis equal; colormap(jet(50));
    scale = 'lin'; lw = 0; title("Umag"); colorbar(); grid off; shading interp;
    fvmplotfield(Umag,scale,lw, 1);
%     saveas(gcf,fullfile(outputpath,'Umag.png'));
    %close all;
    
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