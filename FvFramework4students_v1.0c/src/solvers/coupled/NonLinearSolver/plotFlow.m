function stop = plotFlow(x, casedef)
[p, u, v] = getPUV(casedef,x);

Ux = Field(casedef.dom.allCells,0);     % Velocity [m/s] (x-comp);
set(Ux,u');
Uy = Field(casedef.dom.allCells,0);     % Velocity [m/s] (y-comp);
set(Uy,v');
P =  Field(casedef.dom.allCells,0); 
set(P,p');

figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("Ux"); colorbar();
fvmplotfield(Ux,scale,lw, 1);
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("Uy"); colorbar();
fvmplotfield(Uy,scale,lw, 1);
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1; title("P"); colorbar();
fvmplotfield(P,scale,lw, 1);
close all;
stop=0; % choose stopping criterion here
end