
clc;close all; clear all;
disp('---------------------------------------------------------------------------------------');
disp('Steady state solution of Lid driven cavity flow using Vorticity- Streamfunction Method');
disp('---------------------------------------------------------------------------------------');
%------- Input parameters---------%

Ulid      = 1;                              % Lid Velocity
Lx        = 1;                              % Width
Ly        = 1;                              % Height
Nx        = 41;                             % X- intervals
Ny        = 41;                             % Y- intervals
Re        = 500;                            % Reynolds Number
relaxation= 0.5;                            % Relaxation factor for stability
Niteri    = 5;                              % number of inner iterations
tol       = 0.01;                           % tolerance

%-------Intermediate terms -------%

nu    = (Ulid*Lx)/Re;                       % kinematic viscosity
Dx    = Lx/Nx;
Dy    = Ly/Ny;
Dx2   = 2*Dx;
Dy2   = 2*Dy;
Dxs   = Dx*Dx;
Dys   = Dy*Dy;
beta  = Dxs/Dys;
beta1 = 2*(1+beta);

%--- initialize vorticity and stream function------%
psi     = zeros(Ny+1,Nx+1);                  % Streamfunction
vort    = psi;                               % Vorticity

%---------- Global iterations-------------------------%
err_max=1;
while(err_max>tol)
    save=vort;
%----------------------------------------%
% Jacobi updating of the stream function %
% at the interior nodes Refer Eq (13.1.12)  %
%----------------------------------------%

for iteri=1:Niteri
    
    for i=2:Ny
        for j=2:Nx
            res      = (psi(i,j+1)+psi(i,j-1)+ beta*psi(i+1,j)+beta*psi(i-1,j)+Dxs*vort(i,j))/beta1-psi(i,j);
            psi(i,j) = psi(i,j) + relaxation*res;
        end
    end
    
end

%------------------------------------------------------%
% Compute the vorticity at boundary grid points        %
% using the velocity boundary conditions Refer Eq(13.1.20  )  %
%------------------------------------------------------%

%------ top and bottom walls ---------%

for j=2:Nx
    vort(1,j)    = (7*psi(1,j)-8*psi(2,j)+psi(3,j))/(2*Dys)- 3*Ulid/Dy; % Top Wall
    vort(Ny+1,j) = (7*psi(Ny+1,j)-8*psi(Ny,j) +psi(Ny-1,j))/(2*Dys);      % Bottom Wall
end

%----- left and right walls -----------%

for i=2:Ny
   
    vort(i,1)    = (7*psi(i,1)-8*psi(i,2)+psi(i,3))/(2*Dxs);             % Left Wall
    vort(i,Nx+1) = (7*psi(i,Nx+1)-8*psi(i,Nx)+psi(i,Nx-1))/(2*Dxs);      % Right Wall
    
end

%------------------------------------------%
% compute the velocity at the grid points  %
%         by central differences           %
%------------------------------------------%
ux = zeros(Ny+1,Nx+1);
uy = zeros(Ny+1,Nx+1);
for i=2:Ny
    for j=2:Nx
        ux(i,j) =   (psi(i-1,j)-psi(i+1,j))/Dy2;                 %  dpsi/dy
        uy(i,j) = - (psi(i,j+1)-psi(i,j-1))/Dx2;                 % -dpsi/dx
    end
end

%---------------------------------------------------------------------------------%
% Iterate on Poisson's equation for the vorticity. Refer Eq (13.1.8 and 13.1.20 ) %
%---------------------------------------------------------------------------------%
source = zeros(Ny+1, Nx+1);
for iteri=1:Niteri
    
    for i=2:Ny
        for j=2:Nx
            source(i,j) =  ux(i,j)*(vort(i,j+1)-vort(i,j-1))/Dx2 + uy(i,j)*(vort(i-1,j)-vort(i+1,j))/Dy2;
            source(i,j) = -source(i,j)/nu;
            res = (vort(i,j+1)+vort(i,j-1) + beta*vort(i+1,j)+beta*vort(i-1,j) +Dxs*source(i,j))/beta1-vort(i,j);
            vort(i,j) = vort(i,j) + relaxation*res;
        end
    end
    
end        
%--------Monitor error ---------%
err_max=max(max(abs(vort-save)));
fprintf('err_max = %d\n',err_max);

end 
fprintf('*----CONVERGED-----*\n');
%% Post- Process

x = 0  : Dx : Lx;
y = Ly :-Dy : 0 ;
[xgr,ygr]=meshgrid(x,y);

%---- Get Velocity Field-----%

for i=2:Ny
    for j=2:Nx
        ux(i,j) =   (psi(i-1,j)-psi(i+1,j))/Dy2;   %  dpsi/dy
        uy(i,j) = - (psi(i,j+1)-psi(i,j-1))/Dx2;   % -dpsi/dx
    end
end
V=sqrt(ux.^2+uy.^2);
ux=ux./V;uy=uy./V;


%---- Get Pressure Field by solving  Poisson equation-----%


%---- Streamlines and Velocity Vectors-----%

figure(1)
% Streamlines
contourf(xgr,ygr,psi);colorbar;
caxis([min(min(psi(2:Ny,2:Nx))),max(max(psi(2:Ny,2:Nx)))]);
xlabel('X');ylabel('Y');
title(sprintf('Streamlines + Velocity plot @ Re = %d',Re))
% Velocity vector plot
hold on
quiver(xgr(2:Ny,2:Nx),ygr(2:Ny,2:Nx),ux(2:Ny,2:Nx),uy(2:Ny,2:Nx),0.5,'-k')
hold off



%---- Vorticity contour plot-----%

figure(2)

contourf(xgr,ygr,vort);colorbar;
caxis([min(min(vort(2:Ny,2:Nx))),max(max(vort(2:Ny,2:Nx)))]);
xlabel('X');ylabel('Y');
title(sprintf('Vorticity Contours @ Re = %d',Re))


%----- Centerline Velocity ------%
figure (3)
y_c=zeros(Ny+1,1);
ux_c=y_c;
uy_c=y_c;
for i=1:Ny+1
    y_c(i) = ygr(Ny+2-i,0.5*(Nx+1));
    ux_c(i)= ux (Ny+2-i,0.5*(Nx+1));
    uy_c(i)= uy (Ny+2-i,0.5*(Nx+1));
end
plot(y_c,ux_c,'-k')
hold on
plot(y_c,uy_c,'-g')
xlabel('y-centreline');ylabel('Centreline Velocity');
title(sprintf('Centreline Velocities @ Re = %d',Re));
legend('ux_c','uy_c')
hold off


