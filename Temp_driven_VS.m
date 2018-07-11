
clc;close all; clear all;
disp('---------------------------------------------------------------------------------------');
disp('Steady state solution of Temperature driven cavity flow using Vorticity- Streamfunction Method');
disp('---------------------------------------------------------------------------------------');
%------- Input parameters---------%

Ulid      = 0;                                % Lid Velocity
Lx        = 1;                                % Width
Ly        = 1;                                % Height
Nx        = 31;                               % X- intervals
Ny        = 31;                               % Y- intervals
Pr        =  1;                               % Prandtl Number
Ra        = 1e3;                              % Rayleigh Number
relaxation= 0.2;                              % Relaxation factor for stability
Niteri    = 5;                                % number of inner iterations
tol       = 0.001;                             % tolerance

%-------Intermediate terms -------%
Dx    = Lx/Nx;
Dy    = Ly/Ny;
Dx2   = 2*Dx;
Dy2   = 2*Dy;
Dxs   = Dx*Dx;
Dys   = Dy*Dy;
beta  = Dxs/Dys;
beta1 = 2*(1+beta);

%--- Initialize vorticity and stream function------%
psi     = zeros(Ny+1,Nx+1);                  % Streamfunction + fixed BC for psi
vort    = psi;                               % Vorticity
temp    = psi;                               % Temperature

%% Fixed BC for Temp
temp(2:Ny,1)   =  0;
temp(2:Ny,Nx+1)=  1;


%---------- Global iterations-------------------------%
err_max=1;
while (err_max > tol)

 vort_old = vort;

%--------------------------------------------------%
% Iterate on Poissons equation for stream function %
%--------------------------------------------------%

for iteri=1:Niteri
    
    for i=2:Ny
        for j=2:Nx
            res      = (psi(i,j+1)+psi(i,j-1)+ beta*psi(i+1,j)+beta*psi(i-1,j)+Dxs*vort(i,j))/beta1-psi(i,j);
            psi(i,j) = psi(i,j) + relaxation*res;
        end
    end
    
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
%-----------------------------------------------%
% Iterate on Poisson's equation for temperature %
%-----------------------------------------------%
source_temp = zeros(Ny+1, Nx+1);
for iteri=1:Niteri
    
        for j=2:Nx
            % Top (dtemp/dy=0)
            source_temp(1,j) =  ux(1,j)*(temp(1,j+1)-temp(1,j-1))/Dx2 ;
            source_temp(1,j) = -source_temp(1,j);
            res = (temp(1,j+1)+temp(1,j-1) + 2*beta*temp(2,j)+ Dxs*source_temp(1,j))/beta1-temp(1,j);
            temp(1,j) = temp(1,j) + relaxation*res;
            
            % Bottom (dtemp/dy=0)
            source_temp(Ny+1,j) =  ux(Ny+1,j)*(temp(Ny+1,j+1)-temp(Ny+1,j-1))/Dx2 ;
            source_temp(Ny+1,j) = -source_temp(Ny+1,j);
            res = (temp(Ny+1,j+1)+temp(Ny+1,j-1) + 2*beta*temp(Ny,j)+Dxs*source_temp(Ny+1,j))/beta1-temp(Ny+1,j);
            temp(Ny+1,j) = temp(Ny+1,j) + relaxation*res;
        end
    
    % Sandwich
    for i=2:Ny
        for j=2:Nx
            source_temp(i,j) =  ux(i,j)*(temp(i,j+1)-temp(i,j-1))/Dx2 + uy(i,j)*(temp(i-1,j)-temp(i+1,j))/Dy2;
            source_temp(i,j) = -source_temp(i,j);
            res = (temp(i,j+1)+temp(i,j-1) + beta*temp(i+1,j)+beta*temp(i-1,j) +Dxs*source_temp(i,j))/beta1-temp(i,j);
            temp(i,j) = temp(i,j) + relaxation*res;
        end
    end
  
    
end    

%-------------------------------------------------------%
% Compute the vorticity at boundary grid points         %
% using the velocity boundary conditions. Refer Eq(13.1.20  )  %
%-------------------------------------------------------%

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



%---------------------------------------------------------------%
% Iterate on Poisson's equation for the vorticity. Refer Eq ( ) %
%---------------------------------------------------------------%
source_vort = zeros(Ny+1, Nx+1);
for iteri=1:Niteri
    
    for i=2:Ny
        for j=2:Nx
            source_vort(i,j) =  (ux(i,j)*(vort(i,j+1)-vort(i,j-1))/Dx2 + uy(i,j)*(vort(i-1,j)-vort(i+1,j))/Dy2)/Pr-((temp(i,j+1)-temp(i,j-1))/Dx2)*Ra;
            source_vort(i,j) = -source_vort(i,j);
            res = (vort(i,j+1)+vort(i,j-1) + beta*vort(i+1,j)+beta*vort(i-1,j) +Dxs*source_vort(i,j))/beta1-vort(i,j);
            vort(i,j) = vort(i,j) + relaxation*res;
        end
    end
    
end        
%--------Monitor error ---------%
err_max=max(max(abs(vort-vort_old)));
fprintf('err_max = %d\n',err_max);
end 
fprintf('*--------------------CONVERGED----------------------------*\n');

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


%% Post process

V=sqrt(ux.^2+uy.^2);
ux=ux./V;uy=uy./V;

figure(1)
% Streamlines
contourf(xgr,ygr,psi);colorbar;
caxis([min(min(psi(2:Ny,2:Nx))),max(max(psi(2:Ny,2:Nx)))]);
xlabel('X');ylabel('Y');
title(sprintf('Streamlines + Velocity plot @ Ra = %d and Pr= %d',Ra,Pr))
% Velocity vector plot
hold on
quiver(xgr(2:Ny,2:Nx),ygr(2:Ny,2:Nx),ux(2:Ny,2:Nx),uy(2:Ny,2:Nx),0.5,'-k')
hold off



%---- Vorticity contour plot-----%

figure(2)

contourf(xgr,ygr,vort);colorbar;
caxis([min(min(vort(2:Ny,2:Nx))),max(max(vort(2:Ny,2:Nx)))]);
xlabel('X');ylabel('Y');
title(sprintf('Vorticity Contours @ Ra = %d and Pr = %d',Ra,Pr))

%---- Temperature contour plot-----%
figure(3)

%pcolor(xgr,ygr,temp),shading interp;
contourf(xgr,ygr,temp);
colorbar;
caxis([min(min(temp(2:Ny,2:Nx))),max(max(temp(2:Ny,2:Nx)))]);
xlabel('X');ylabel('Y');
title(sprintf('Temperature Contours @ Ra = %d and Pr = %d',Ra, Pr))

%----- Centerline Velocity ------%
figure (4)
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
title(sprintf('Centreline Velocities @ Ra = %d and Pr= %d',Ra, Pr));
legend('ux_c','uy_c')
hold off



