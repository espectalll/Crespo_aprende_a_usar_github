% TIERRA HUECA: Campo Gravitacional de una Corteza esférica
% Hecho deprisa y corriendo (i.e. muy cutre) por QuantumFracture
 
clc; clear all; close all;
 
XX=2; %PRECISION DEL CALCULO
 
rf=6371e3; %Radio Externo (m)
r0=0.8*rf; %Radio Interno
omf=pi/2*0.2; %Tamaño angular del agujero
 
L=r0; %Dimensión
h=L/6; %Pasos
n=length(-L:h:L);
hr=(rf-r0)/3/XX;
hphi=0.3/XX;
hom=0.3/XX;
 
G=6.67e-11; %N m^2/kg^2
 
% Grid en x,y,z
[x, y, z] = meshgrid(-L:h:L, -L:h:L, -L:h:L);
 
% Posiciones de la Corteza
 
R=[];
num=0;
 
for r=r0:hr:rf
    for phi=0:hphi:(2*pi)
        for om=omf:hom:(pi-omf)
            num=num+1;
            x0=r*sin(om)*cos(phi);
            y0=r*sin(om)*sin(phi);
            z0=r*cos(om);
            hold on
            plot3(x0,y0,z0,'ok') %Plot de la Corteza
            R=[R [x0; y0; z0;] ];
        end
    end
end
 
M=5.97e24*(rf^3-r0^3)/rf^3/num; %kg %Masa de cada unidad de corteza
 
F=[0;0;0];
 
for i=1:n
    for j=1:n
        for k=1:n
            if norm([x(i,j,k); y(i,j,k);z(i,j,k)])>r0
                F=[0;0;0];
            else
                for q=1:length(R)
                vec=[x(i,j,k); y(i,j,k);z(i,j,k)] - R(:,q);
                vec=vec/norm(vec);
                f=-G*M/(norm([x(i,j,k); y(i,j,k);z(i,j,k)] - R(:,q)))^2*vec;
                F=F+f;
            end
            
            u(i,j,k)=F(1);
            v(i,j,k)=F(2);
            w(i,j,k)=F(3);
            F=0;
            end
        end
    end
end
 
%Máxima Gravedad en el Borde
maxim1=mean(mean(Mod(:,:,2)));
maxim2=mean(mean(Mod(:,:,end-1)));
grav=max([maxim1 maxim2])
 
% Plot del Campo Gravitacional
quiver3(x, y, z, u, v, w)
hold off
% Limites (ajustados al interior de la corteza)
axis([-L L -L L -L L])
 
% Titulo y labels
title('Gravedad Corteza Esférica')
xlabel('x')
ylabel('y')
zlabel('z')
