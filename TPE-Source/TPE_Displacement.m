function [U,V,W] = TPE_Displacement(x_obs,y_obs,z_obs,tpe_x,tpe_y,c,a,d,Delta_p,Delta_T,H,alpha,nu,mu,varargin)

% TPE_Displacement: Thermo-Poro-Elastic (Source) Displacement. Computes the 
%               displaceent components at points (x_obs, y_obs, z_obs) in a 
%               half-space, due to a TPE centred at point (tpe_x, tpe_y, c).          
% 
% Arguments: (input)
% 
%      x_obs, y_obs, z_obs  - The observation point location in cartesian coordinates. 
%      They should be provided either as a NxN grid, or as a
%      one-dimensional array (e.g. x_obs=[x_obs(1)...x_obs(N)],y_obs=0)
% 
%      tpe_x, tpe_y - The horizontal cartesian coordinates of the TPE
%      centre.
% 
%      c - The depth of the TPE centre with respect to the free surface
%      (must be positive).
% 
%      a - The radius of the TPE
% 
%      d - The thickness of the TPE
% 
%      Delta_p - The increase in pore pressure within the TPE (must be
%      positive)
% 
%      Delta_T - The increase in temperature within the TPE (must be
%      positive)
% 
%      H - Biot's constant for the medium 
%      (see Mantiloni et al., 2020, Equation 3)
%         
%      alpha - coefficient of thermal expansion for the medium 
%      (see Mantiloni et al., 2020, Equation 3)
%         
%      nu - Poisson's ratio for the medium 
%         
%      mu - Modulus of rigidity for the medium 
% 
%      varargin: % n - The number of coefficients in the Legendre's 
%      polynomials expansion for singular displacement components (default is 
%      200), see Mantiloni et al., 2020, Equation 12.
%     
% Arguments: (output)
% 
%      U,V,W - The displacement components in a cartesian reference frame at point 
%      (x_obs,y_obs,z_obs). U is the x-component, V the y-component, W the
%      vertical component. For NxN obs points, there will be three Nx1
%      arrays.
% 
% Example usage:
% 
%  x_obs = linspace(-1,1,100)*1e4;
%  y_obs = 0;
%  z_obs = 0;
%  tpe_x = 0; tpe_y = 0;
%  c = 3e3;
%  a = 1e3;
%  d = 50;
%  Delta_p = 1e3;
%  Delta_T = 100;
%  H = 1e11;
%  alpha = 3e-5;
%  nu = 0.3;
%  mu = 6e9;
%  n = 200;
%  [U,V,W] = TPE_Displacement(x_obs,y_obs,z_obs,tpe_x,tpe_y,c,a,d,Delta_p,Delta_T,H,alpha,nu,mu,n);
%  figure
%  ax=axes; 
%  hold on;
%  plot(x_obs,U,'LineWidth',1.5)
%  plot(x_obs,V,'LineWidth',1.5)
%  plot(x_obs,W,'LineWidth',1.5);
%  SP1=-a; 
%  line([SP1 SP1],get(ax,'YLim'),'Color','k','LineWidth',2)
%  SP2=a; 
%  line([SP2 SP2],get(ax,'YLim'),'Color','k','LineWidth',2)
%  legend('U','V','W','TPE boundary');
%  xlabel('x (m)')
%  ylabel('Displacement (m)')
%  hold off
% 
% Author: Lorenzo Mantiloni
% University of Potsdam / German Research Centre for Geosciences GFZ / Università di Bologna
% 
% Copyright (C) 2019  Lorenzo Mantiloni
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%Checking inputs
if ~isempty(varargin)
    n = varargin{1,1};
else
    n = 200;
end

if d/a > 0.1
    error('Warning: the TPE aspect ratio is too high and solutions will fail')
%     return
end 
if Delta_p < 0
    error('Pore pressure increase cannot be negative')
%     return
end 
if Delta_T < 0
    error('Temperature increase cannot be negative')
%     return
end 

if size(x_obs,2) > 1
    x_obs = x_obs(:);
end
if size(y_obs,2) > 1
    y_obs = y_obs(:);
end
if size(x_obs,1) ~= size(y_obs,1)
    if numel(x_obs) == 1 
        x_obs = y_obs*0 + x_obs;
    else
        if numel(y_obs) == 1 
            y_obs = x_obs*0 + y_obs;
        else
            error('Please provide observation points either as a NxN grid or as [x_obs(1)...x_obs(N)],y_obs or viceversa'\n)
        end
    end
end
        
if numel(z_obs) > 1
    error('Please select a common depth for all observation points'\n)
end
p = numel(x_obs);    
c1=c-z_obs; %Vertical distance between TPE median plane and observation point

%Pre-allocating some arrays
U=zeros(p,1);
V=zeros(p,1);
W=zeros(p,1);

%Setting some parameters
lambda=2*mu*nu/(1-2*nu); %Lamé's first parameter of medium
k=(lambda + (2/3)*mu); %Bulk modulus of medium
e0 = (1/(3*H))*Delta_p + (1/3)*alpha*Delta_T;
ke0=k*e0;

for count=1:p
    x_o = x_obs(count);
    y_o = y_obs(count);

%Preparing the cartesian-to-spherical coordinate transformation for the
%solution of singular components

R = ((x_o-tpe_x).^2 + (y_o-tpe_y).^2).^(1/2); %Radial coordinates of observation points

%Azimuthal angle of observation points
if c1~=0
   TH = atan(R./c1);
else
   TH = pi/2;
end

%Polar angle of observation points
if (x_o-tpe_x)<0 
   PHI = atan((y_o-tpe_y)./(x_o-tpe_x)) + pi;
else
   PHI = atan((y_o-tpe_y)./(x_o-tpe_x));
end

%Initializing some arrays
coeff=zeros(2*n,1);

%Setting some constants
C=1/(16*pi*mu*(1-nu));
e1=e0*(1+nu)/(1-nu);
A=e1*(d/(2*a));
K1=2*C*(1-2*nu)*(3-4*nu)*3*ke0;
K2=4*C*(1-2*nu)*3*ke0;

%This is to avoid oddities in TH and PHI
ctnan = isnan(TH);
TH(ctnan) = 0;
cpnan = isnan(PHI);
PHI(cpnan) = 0;

%Retrieving the coefficients of Legendre's polynomials
coeff(1)=1;
for i=3:2*n
    coeff(i)=-((i-2)/(i-1))*coeff(i-2);
end

%Calculating the series' terms and the singular components in spherical
%coordinates
    Pol=Leg_Pol(TH,2*n);
    dPol=Leg_Pol_Der(TH,2*n,Pol);
    x=cos(TH);
    r=(R^2 + c1^2)^(1/2); 
    P=Pol_Leg_Sum(r,a,n,coeff,Pol);
    dP=Pol_Leg_Derivative_Sum(r,a,n,coeff,dPol);
        if r<=a       
            u_r= A*a*(abs(x) + P);
            u_t= -A*a*(sign(x)*sin(TH) - dP);
        else
            u_r= A*a*P;
            u_t= -A*a*dP;
        end    

%Converting back to cartesian coordinates
        u_x= +sin(TH)*cos(PHI)*u_r + cos(TH)*cos(PHI)*u_t;
        u_y= +sin(TH)*sin(PHI)*u_r + cos(TH)*sin(PHI)*u_t;
        u_z= -cos(TH)*u_r + sin(TH)*u_t; %Bear in mind that the z-axis is positive downward

%Retrieving the non-singular components and the complete ones

%First component
    fun_a = @(y) log((sqrt((x_o-sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2+(y_o-y).^2+(z_obs+c+d/2).^2) + z_obs + c + d/2)./(sqrt((x_o-sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2+(y_o-y).^2+(z_obs+c-d/2).^2) + z_obs + c - d/2)) - log((sqrt((x_o+sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2+(y_o-y).^2+(z_obs+c+d/2).^2) + z_obs + c + d/2)./(sqrt((x_o+sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2+(y_o-y).^2+(z_obs+c-d/2).^2) + z_obs + c - d/2));
    fun_b = @(y) 1./sqrt((x_o-sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2 +(y_o-y).^2+(z_obs+c+d/2).^2) - 1./sqrt((x_o+sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2 +(y_o-y).^2+(z_obs+c+d/2).^2) - 1./sqrt((x_o-sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2 +(y_o-y).^2+(z_obs+c-d/2).^2) + 1./sqrt((x_o+sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2 +(y_o-y).^2+(z_obs+c-d/2).^2);
    q_a = integral(fun_a,-a+tpe_y,a+tpe_y,'RelTol',0,'AbsTol',1e-6);
    q_b = integral(fun_b,-a+tpe_y,a+tpe_y,'RelTol',0,'AbsTol',1e-6);
    int_single1= K1*q_a + K2*z_obs.*q_b + u_x;  

%Second component
        fun_a = @(x) log((sqrt((y_o-sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2+(x_o-x).^2+(z_obs+c+d/2).^2) + z_obs + c + d/2)./(sqrt((y_o-sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2+(x_o-x).^2+(z_obs+c-d/2).^2) + z_obs + c - d/2)) - log((sqrt((y_o+sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2+(x_o-x).^2+(z_obs+c+d/2).^2) + z_obs + c + d/2)./(sqrt((y_o+sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2+(x_o-x).^2+(z_obs+c-d/2).^2) + z_obs + c - d/2));
        fun_b = @(x) 1./sqrt((y_o-sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2 +(x_o-x).^2+(z_obs+c+d/2).^2) - 1./sqrt((y_o+sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2 +(x_o-x).^2+(z_obs+c+d/2).^2) - 1./sqrt((y_o-sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2 +(x_o-x).^2+(z_obs+c-d/2).^2) + 1./sqrt((y_o+sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2 +(x_o-x).^2+(z_obs+c-d/2).^2);
        q_a = integral(fun_a,-a+tpe_x,a+tpe_x,'RelTol',0,'AbsTol',1e-6);
        q_b = integral(fun_b,-a+tpe_x,a+tpe_x,'RelTol',0,'AbsTol',1e-6);
        int_single2= K1*q_a + K2*z_obs.*q_b + u_y;

%Third component
        fun_a = @(y) log((sqrt((x_o-sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2+(y_o-y).^2+(z_obs+c+d/2).^2) - x_o + sqrt(a^2 - (y-tpe_y).^2)+tpe_x)./(sqrt((x_o+sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2+(y_o-y).^2+(z_obs+c+d/2).^2) - x_o - sqrt(a^2 - (y-tpe_y).^2)+tpe_x)) - log((sqrt((x_o-sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2+(y_o-y).^2+(z_obs+c-d/2).^2) - x_o + sqrt(a^2 - (y-tpe_y).^2)+tpe_x)./(sqrt((x_o+sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2+(y_o-y).^2+(z_obs+c-d/2).^2) - x_o - sqrt(a^2 - (y-tpe_y).^2)+tpe_x));
        fun_b = @(y) (-(z_obs + c + d/2).*((x_o - sqrt(a^2 - (y-tpe_y).^2)+tpe_x)./(sqrt((x_o-sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2+(y_o-y).^2+(z_obs+c+d/2).^2).*((y_o - y).^2 + (z_obs + c + d/2).^2)) - (x_o + sqrt(a^2 - (y-tpe_y).^2)+tpe_x)./(sqrt((x_o+sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2+(y_o-y).^2+(z_obs+c+d/2).^2).*((y_o - y).^2 + (z_obs + c + d/2).^2))) + (z_obs + c - d/2).*((x_o - sqrt(a^2 - (y-tpe_y).^2)+tpe_x)./(sqrt((x_o-sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2+(y_o-y).^2+(z_obs+c-d/2).^2).*((y_o - y).^2 + (z_obs + c - d/2).^2)) - (x_o + sqrt(a^2 - (y-tpe_y).^2)+tpe_x)./(sqrt((x_o+sqrt(a^2 - (y-tpe_y).^2)-tpe_x).^2+(y_o-y).^2+(z_obs+c-d/2).^2).*((y_o - y).^2 + (z_obs + c - d/2).^2))));
        q_a = integral(fun_a,-a+tpe_y,a+tpe_y,'RelTol',0,'AbsTol',1e-6);
        q_b = integral(fun_b,-a+tpe_y,a+tpe_y,'RelTol',0,'AbsTol',1e-6);
        int_single3= K1*q_a + K2*z_obs.*q_b + u_z;    

u=int_single1;
v=int_single2;
w=int_single3;

U(count)=u;
V(count)=v;
W(count)=-w; %The vertical axis is positive downward in Mantiloni et al., 2020 

end

end

%Defining some functions retrieving the Legendre's polynomials
function Pol=Leg_Pol(v,l)
Pol=zeros(l,1);
Pol(1)=1;
Pol(2)=cos(v);
if l>1
    for i=2:l-1
        Pol(i+1)=((2*(i-1)+1)/i)*cos(v)*Pol(i) - ((i-1)/i)*Pol(i-1);
    end
end
end

function dPol=Leg_Pol_Der(v,q,Pol)
dPol_aux=zeros(q,1);
dPol_aux(1)=0;
dPol_aux(2)=1;
if q>1
    for i=2:q-1
        dPol_aux(i+1)=(2*(i-1)+1)*Pol(i) + dPol_aux(i-1);
    end
end
dPol=-sin(v)*dPol_aux;
end

function P=Pol_Leg_Sum(r,a,n,coeff,Pol)
P=0;
    if r<=a
        for j=2:n
            Q=P;
            P=Q + coeff(2*j-1)*Pol(2*j-1)*(2*(j-1)/(2*(j-1)-1))*(r/a)^(2*(j-1)-1);
        end
    else
        for j=1:n
            Q=P;
            P=Q + coeff(2*j-1)*Pol(2*j-1)*((a/r)^(2*(j-1)+2))*((2*(j-1)+1)/(2*(j-1)+2));
        end
    end
end

function dP = Pol_Leg_Derivative_Sum(r,a,n,coeff,dPol)
dP=0;
    if r<=a
        for j=2:n
            Q=dP;
            dP=Q + coeff(2*j-1)*dPol(2*j-1)*(1/(2*(j-1)-1))*((r/a)^(2*(j-1)-1));
        end
    else
        for j=2:n
            Q=dP;
            dP=Q + coeff(2*j-1)*dPol(2*j-1)*((a/r)^(2*(j-1)+2))*(1/(2*(j-1)+2));
        end
    end
end
