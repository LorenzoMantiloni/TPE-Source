function [Stress,Strain] = TPE_Stresses(x_obs,y_obs,z_obs,tpe_x,tpe_y,c,a,d,Delta_p,Delta_T,H,alpha,nu,mu,varargin)

% TPE_Stresses: Thermo-Poro-Elastic (Source) Stresses. Computes the 
%               components of the stress tensor at points (x_obs, y_obs,
%               z_obs) in a half-space, due to a TPE centred at point
%               (tpe_x, tpe_y, c).
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
%      Stress - The stress tensor in a cartesian reference frame at point 
%      (x_obs,y_obs,z_obs), given as a one-row array for each obs point
%      [T_xx,T_yy,T_zz,T_xy,T_xz,T_yz] (for N obs points, it will populate a
%      N x 6 matrix)
% 
%      Strain - The strain tensor in a cartesian reference frame at point 
%      (x_obs,y_obs,z_obs) given as a one-row array for each os point
%      [E_xx,E_yy,E_zz,E_xy,E_xz,E_yz] (for N obs points, it will populate a
%      N x 6 matrix)
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
%  [Stress,Strain] = TPE_Stresses(x_obs,y_obs,z_obs,tpe_x,tpe_y,c,a,d,Delta_p,Delta_T,H,alpha,nu,mu,n);
%  figure
%  ax=axes; 
%  hold on;
%  plot(x_obs,Stress(:,1),'LineWidth',1.5)
%  plot(x_obs,Stress(:,2),'LineWidth',1.5)
%  plot(x_obs,Stress(:,3),'LineWidth',1.5);
%  plot(x_obs,Stress(:,5),'LineWidth',1.5);
%  SP1=-a; 
%  line([SP1 SP1],get(ax,'YLim'),'Color','k','LineWidth',2)
%  SP2=a; 
%  line([SP2 SP2],get(ax,'YLim'),'Color','k','LineWidth',2)
%  legend('\sigma_{xx}','\sigma_{yy}','\sigma_{zz}','\sigma_{xz}','TPE boundary');
%  xlabel('x (m)')
%  ylabel('\sigma (Pa)')
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

%Pre-allocating output arrays
Stress=zeros(p,6);
Strain=zeros(p,6);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Strain tensor components:

%Calculating the series' terms and the singular components in spherical
%coordinates
    Pol=Leg_Pol(TH,2*n);
    dPol=Leg_Pol_Der(TH,2*n,Pol);
    ddPol=Leg_Pol_2Der(TH,2*n,Pol,dPol);
    r=(R^2 + c1^2)^(1/2); 
    P1=Pol_Leg_Sum(r,a,n,1,coeff,Pol);
    P2=Pol_Leg_Sum(r,a,n,2,coeff,Pol);
    ddP2=Pol_Leg_2Derivative_Sum(r,a,n,2,coeff,ddPol);
    P3=Pol_Leg_Sum(r,a,n,3,coeff,Pol);
    dP3=Pol_Leg_Derivative_Sum(r,a,n,3,TH,coeff,dPol);
    dP4=Pol_Leg_Derivative_Sum(r,a,n,4,TH,coeff,dPol);
    if r<=a       
        e_rr= A*P1;
        e_tt= A*(2*dirac(cos(TH))*(a/r) + P2 + ddP2);
        e_ff= A*(dP3 + P3);
        e_rt= A*dP4;
        %Converting back to cartesian coordinates
        e_zz= e_rr*(cos(TH))^2 - 2*e_rt*sin(TH)*cos(TH) + e_tt(i)*(sin(TH))^2;
        e_xz= e_rr*cos(TH)*sin(TH)*cos(PHI) + e_rt*cos(PHI)*((cos(TH))^2 - (sin(TH))^2) - e_tt*cos(TH)*sin(TH)*cos(PHI);
        e_xx= e_rr*((sin(TH))^2)*(cos(PHI)^2) + 2*e_rt*sin(TH)*cos(TH)*(cos(PHI))^2 + e_tt*((cos(TH))^2)*(cos(PHI))^2 + e_ff*(sin(PHI))^2;
        e_yy= e_rr*((sin(PHI))^2)*(sin(TH))^2 + 2*e_rt*cos(TH)*sin(TH)*(sin(PHI))^2 + e_tt*((cos(TH))^2)*(sin(PHI))^2 + e_ff*(cos(PHI))^2;
    else
        e_rr= -A*P1;
        e_tt= (A/2)*(((a/r)^3) + P2 - ddP2 );
        e_ff= -A*(dP3 - P3);
        e_rt= A*dP4;
        %Converting back to cartesian coordinates
        e_zz= e_rr*(cos(TH))^2 - 2*e_rt*sin(TH)*cos(TH) + e_tt*(sin(TH))^2;
        e_xz= e_rr*cos(TH)*sin(TH)*cos(PHI) + e_rt*cos(PHI)*((cos(TH))^2 - (sin(TH))^2) - e_tt*cos(TH)*sin(TH)*cos(PHI);
        e_xx= e_rr*((sin(TH))^2)*(cos(PHI)^2) + 2*e_rt*sin(TH)*cos(TH)*(cos(PHI))^2 + e_tt*((cos(TH))^2)*(cos(PHI))^2 + e_ff*(sin(PHI))^2;
        e_yy= e_rr*((sin(PHI))^2)*(sin(TH))^2 + 2*e_rt*cos(TH)*sin(TH)*(sin(PHI))^2 + e_tt*((cos(TH))^2)*(sin(PHI))^2 + e_ff*(cos(PHI))^2;
    end
    if isnan(e_xx)
       e_xx = 0;
    end
    if isnan(e_yy)
       e_yy = 0;
    end
    if isnan(e_zz)
       e_zz = 0;
    end
    if ((c - d/2) < z_obs) && (z_obs < (c + d/2)) && (abs(x_o - tpe_x) < a) && (abs(y_o - tpe_y) < a)
        sum_e = e_zz + e_xx + e_yy + e1;
    else 
        sum_e = e_zz + e_xx + e_yy;
    end

%Singular stresses:

    if ((c-d/2) < z_obs) && (z_obs < c+d/2) && (x_o < a)
        ts_zz= 2*mu*e_zz;%/(2*mu*e1);
        ts_xx= 2*mu*(-e1 + e_xx);%/(2*mu*e1);
        ts_yy= 2*mu*(-e1 + e_yy);%/(2*mu*e1);
        tsd_xz= -2*mu*e_xz;
    else
        ts_zz= (2*mu*e_zz + lambda*sum_e);%/(2*mu*e1);
        ts_xx= (2*mu*e_xx + lambda*sum_e);%/(2*mu*e1);
        ts_yy= (2*mu*e_yy + lambda*sum_e);%/(2*mu*e1);
        tsd_xz= -2*mu*e_xz;
    end

% % t_xy:
%     fun_a = @(z) 1./sqrt(x_o^2 +(y_o-a).^2+(z_obs+z).^2) - 1./sqrt(x_o^2 +(y_o+a).^2+(z_obs+z).^2) - 1./sqrt(x_o^2 +(y_o-a).^2+(z_obs+z).^2) + 1./sqrt(x_o^2 +(y_o+a).^2+(z_obs+z).^2);
%     fun_b = @(z) 1./(sqrt(x_o^2 +(y_o-a).^2+(z_obs+z).^2)).^3 - 1./(sqrt(x_o^2 +(y_o+a).^2+(z_obs+z).^2)).^3 - 1./(sqrt(x_o^2 +(y_o-a).^2+(z_obs+z).^2)).^3 + 1./(sqrt(x_o^2 +(y_o+a).^2+(z_obs+z).^2)).^3;
%     q_a = integral(fun_a,c-d/2,c+d/2,'RelTol',1e-12,'AbsTol',1e-12);
%     q_b = integral(fun_b,c-d/2,c+d/2,'RelTol',1e-12,'AbsTol',1e-12);
%     T_xy= K1*q_a - K2*z_obs.*q_b; 



%Partial derivatives of displacement components: 
%du_x/dx
    fun_a = @(y) (x_o - sqrt(a^2 - (y-tpe_y).^2) - tpe_x).*(((z_obs + c + d/2)./(((y_o - y).^2 + (x_o- sqrt(a^2 - (y-tpe_y).^2) - tpe_x).^2).*(sqrt((x_o-sqrt(a^2 - (y-tpe_y).^2) - tpe_x).^2+(y_o-y).^2+(z_obs+c+d/2).^2)))) - ((z_obs + c - d/2)./(((y_o - y).^2 + (x_o- sqrt(a^2 - (y-tpe_y).^2) - tpe_x).^2).*(sqrt((x_o-sqrt(a^2 - (y-tpe_y).^2) - tpe_x).^2+(y_o-y).^2+(z_obs+c-d/2).^2))))) - (x_o + sqrt(a^2 - (y-tpe_y).^2) - tpe_x).*(((z_obs + c + d/2)./(((y_o - y).^2 + (x_o+ sqrt(a^2 - (y-tpe_y).^2) - tpe_x).^2).*(sqrt((x_o+sqrt(a^2 - (y-tpe_y).^2) - tpe_x).^2+(y_o-y).^2+(z_obs+c+d/2).^2)))) - ((z_obs + c - d/2)./(((y_o - y).^2 + (x_o+ sqrt(a^2 - (y-tpe_y).^2) - tpe_x).^2).*(sqrt((x_o+sqrt(a^2 - (y-tpe_y).^2) - tpe_x).^2+(y_o-y).^2+(z_obs+c-d/2).^2))))); 
    fun_b = @(y) (x_o - sqrt(a^2 - (y-tpe_y).^2) - tpe_x)./(sqrt((x_o-sqrt(a^2 - (y-tpe_y).^2) - tpe_x).^2+(y_o-y).^2+(z_obs+c+d/2).^2).^3) - (x_o - sqrt(a^2 - (y-tpe_y).^2) - tpe_x)./(sqrt((x_o-sqrt(a^2 - (y-tpe_y).^2) - tpe_x).^2+(y_o-y).^2+(z_obs+c-d/2).^2).^3) - (x_o + sqrt(a^2 - (y-tpe_y).^2) - tpe_x)./(sqrt((x_o+sqrt(a^2 - (y-tpe_y).^2) - tpe_x).^2+(y_o-y).^2+(z_obs+c+d/2).^2).^3) + (x_o + sqrt(a^2 - (y-tpe_y).^2) - tpe_x)./(sqrt((x_o+sqrt(a^2 - (y-tpe_y).^2) - tpe_x).^2+(y_o-y).^2+(z_obs+c-d/2).^2).^3);
    q_a = integral(fun_a,-a+tpe_y,a+tpe_y,'RelTol',1e-12,'AbsTol',1e-12);
    q_b = integral(fun_b,-a+tpe_y,a+tpe_y,'RelTol',1e-12,'AbsTol',1e-12);
    du_x= -K1*q_a - K2*z_obs.*q_b;

%du_y/dy
    fun_a = @(x) (y_o - sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).*((z_obs + c + d/2)./(((x_o - x).^2 + (y_o- sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2).*(sqrt((y_o-sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2+(x_o-x).^2+(z_obs+c+d/2).^2))) - (z_obs + c - d/2)./(((x_o - x).^2 + (y_o- sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2).*(sqrt((y_o-sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2+(x_o-x).^2+(z_obs+c-d/2).^2)))) - (y_o + sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).*((z_obs + c + d/2)./(((x_o - x).^2 + (y_o+ sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2).*(sqrt((y_o+sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2+(x_o-x).^2+(z_obs+c+d/2).^2))) - (z_obs + c - d/2)./(((x_o - x).^2 + (y_o+ sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2).*(sqrt((y_o+sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2+(x_o-x).^2+(z_obs+c-d/2).^2)))); 
    fun_b = @(x) (y_o - sqrt((a)^2 - (x-tpe_x).^2)-tpe_y)./(sqrt((y_o-sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2+(x_o-x).^2+(z_obs+c+d/2).^2).^3) - (y_o - sqrt((a)^2 - (x-tpe_x).^2)-tpe_y)./(sqrt((y_o-sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2+(x_o-x).^2+(z_obs+c-d/2).^2).^3) - (y_o + sqrt((a)^2 - (x-tpe_x).^2)-tpe_y)./(sqrt((y_o+sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2+(x_o-x).^2+(z_obs+c+d/2).^2).^3) + (y_o + sqrt((a)^2 - (x-tpe_x).^2)-tpe_y)./(sqrt((y_o+sqrt((a)^2 - (x-tpe_x).^2)-tpe_y).^2+(x_o-x).^2+(z_obs+c-d/2).^2).^3);
    q_a = integral(fun_a,-a+tpe_x,a+tpe_x,'RelTol',1e-12,'AbsTol',1e-12);
    q_b = integral(fun_b,-a+tpe_x,a+tpe_x,'RelTol',1e-12,'AbsTol',1e-12);
    du_y= -K1*q_a - K2*z_obs.*q_b;

%du_z/dz 
    fun_a = @(y) -(z_obs + c + d/2).*(-(x_o - (a^2 - (y-tpe_y).^2).^(1/2)-tpe_x)./(((y_o-y).^2 + (z_obs + c + d/2).^2).*((x_o-(a^2 - (y-tpe_y).^2).^(1/2)-tpe_x).^2+(y_o-y).^2+(z_obs+c+d/2).^2).^(1/2)) + (x_o + (a^2 - (y-tpe_y).^2).^(1/2)-tpe_x)./(((y_o-y).^2 + (z_obs + c + d/2).^2).*((x_o+(a^2 - (y-tpe_y).^2).^(1/2)-tpe_x).^2+(y_o-y).^2+(z_obs+c+d/2).^2).^(1/2))) - (z_obs + c - d/2).*((x_o - (a^2 - (y-tpe_y).^2).^(1/2)-tpe_x)./(((y_o-y).^2 + (z_obs + c - d/2).^2).*((x_o-(a^2 - (y-tpe_y).^2).^(1/2)-tpe_x).^2+(y_o-y).^2+(z_obs+c-d/2).^2).^(1/2)) - (x_o + (a^2 - (y-tpe_y).^2).^(1/2)-tpe_x)./(((y_o-y).^2 + (z_obs + c - d/2).^2).*((x_o+(a^2 - (y-tpe_y).^2).^(1/2)-tpe_x).^2+(y_o-y).^2+(z_obs+c-d/2).^2).^(1/2)));
    fun2_a= @(y,x) (z_obs+c+d/2)./(((x_o-x).^2 + (y_o-y).^2 + (z_obs+c+d/2).^2).^(3/2)) - (z_obs+c-d/2)./(((x_o-x).^2 + (y_o-y).^2 + (z_obs+c-d/2).^2).^(3/2)) + z_obs./(((x_o-x).^2 + (y_o-y).^2 + (z_obs+c+d/2).^2).^(3/2)) - z_obs./(((x_o-x).^2 + (y_o-y).^2 + (z_obs+c-d/2).^2).^(3/2)) - 3*z_obs.*((z_obs+c+d/2).^2)./(((x_o-x).^2 + (y_o-y).^2 + (z_obs+c+d/2).^2).^(5/2)) + 3*z_obs.*((z_obs+c-d/2).^2)./(((x_o-x).^2 + (y_o-y).^2 + (z_obs+c-d/2).^2).^(5/2));
    xmin = @(y) -sqrt(a^2 - (y-tpe_y).^2)+tpe_x;
    xmax = @(y) sqrt(a^2 - (y-tpe_y).^2)+tpe_x;
    q_a = integral(fun_a,-a+tpe_y,a+tpe_y,'RelTol',1e-12,'AbsTol',1e-12);
    q2_b= integral2(fun2_a,-a+tpe_y,a+tpe_y,xmin,xmax,'RelTol',1e-12,'AbsTol',1e-12);
    du_z= K1*q_a + K2*q2_b;
    
%e_xzns:
    fun_b = @(y) -z_obs*(((z_obs+c+(d/2))./(((x_o-sqrt(a^2 - y.^2)).^2 + (y_o-y).^2 + (z_obs+c+d/2).^2).^(3/2))) - ((z_obs+c-(d/2))./(((x_o-sqrt(a^2 - y.^2)).^2 + (y_o-y).^2 + (z_obs+c-d/2).^2).^(3/2))) - ((z_obs+c+(d/2))./(((x_o+sqrt(a^2 - y.^2)).^2 + (y_o-y).^2 + (z_obs+c+d/2).^2).^(3/2))) + ((z_obs+c-(d/2))./(((x_o+sqrt(a^2 - y.^2)).^2 + (y_o-y).^2 + (z_obs+c-d/2).^2).^(3/2)))); 
    fun_a = @(y) 1./(((x_o-sqrt(a^2 - y.^2)).^2 + (y_o-y).^2 + (z_obs+c+d/2).^2).^(1/2)) - 1./(((x_o-sqrt(a^2 - y.^2)).^2 + (y_o-y).^2 + (z_obs+c-d/2).^2).^(1/2)) - 1./(((x_o+sqrt(a^2 - y.^2)).^2 + (y_o-y).^2 + (z_obs+c+d/2).^2).^(1/2)) + 1./(((x_o+sqrt(a^2 - y.^2)).^2 + (y_o-y).^2 + (z_obs+c-d/2).^2).^(1/2)) - z_obs*(((z_obs+c+(d/2))./(((x_o-sqrt(a^2 - y.^2)).^2 + (y_o-y).^2 + (z_obs+c+d/2).^2).^(3/2))) - ((z_obs+c-(d/2))./(((x_o-sqrt(a^2 - y.^2)).^2 + (y_o-y).^2 + (z_obs+c-d/2).^2).^(3/2))) - ((z_obs+c+(d/2))./(((x_o+sqrt(a^2 - y.^2)).^2 + (y_o-y).^2 + (z_obs+c+d/2).^2).^(3/2))) + ((z_obs+c-(d/2))./(((x_o+sqrt(a^2 - y.^2)).^2 + (y_o-y).^2 + (z_obs+c-d/2).^2).^(3/2))));  
    q_a = integral(fun_a,-a,a,'RelTol',1e-12,'AbsTol',1e-12);
    q_b = integral(fun_b,-a,a,'RelTol',1e-12,'AbsTol',1e-12);
    e_xzns= (K2/2)*(q_a + q_b);
    
%Divergence of displacement field:
sum_u = du_x + du_y + du_z;

%Non-singular diagonal stresses:
T_xx = (lambda).*sum_u + 2*mu*du_x;
T_yy = (lambda).*sum_u + 2*mu*du_y;
T_zz = (lambda).*sum_u + 2*mu*du_z;

%Non-singular deviatoric stresses:
T_xy=0;
T_xz=2*mu*e_xzns;
T_yz=0;

%Total strain:
E_xx= e_xx + du_x;
E_yy= e_yy + du_y;
E_zz= e_zz + du_z;
E_xy= 0;
E_xz= e_xz + e_xzns;
E_yz= 0;

%Total stresses:
T_xx= T_xx + ts_xx;
T_yy= T_yy + ts_yy;
T_zz= T_zz + ts_zz;
T_xz = T_xz + tsd_xz;

Stress(count,:) = [T_xx,T_yy,T_zz,T_xy,T_xz,T_yz];

Strain(count,:) = [E_xx,E_yy,E_zz,E_xy,E_xz,E_yz];
        
end

end


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

function dPol=Leg_Pol_Der(v,l,Pol)
dPol_aux=zeros(l,1);
dPol_aux(1)=0;
dPol_aux(2)=1;
if l>1
    for i=2:l-1
        dPol_aux(i+1)=(2*(i-1)+1)*Pol(i) + dPol_aux(i-1);
    end
end
dPol=-sin(v)*dPol_aux;
end

function ddPol=Leg_Pol_2Der(v,l,Pol,dPol)
ddPol=zeros(l,1);
ddPol_aux=zeros(l,1);
for i=1:l
    ddPol_aux(i)= 2*cos(v)/((sin(v))^2)*(dPol(i)./(-sin(v))) - ((i-1)*i/((sin(v))^2))*Pol(i);
    ddPol(i)=(sin(v))^2*ddPol_aux(i) + (cos(v)/sin(v))*dPol(i);
end
end

function P=Pol_Leg_Sum(r,a,n,var,coeff,Pol)
P=0;
if var==1
    if r<=a
        for j=2:n
            Q=P;
            P=Q + coeff(2*j-1)*Pol(2*j-1)*((r/a)^(2*(j-1)-2))*2*(j-1);
        end
    else
        for j=1:n
            Q=P;
            P=Q + coeff(2*j-1)*Pol(2*j-1)*((a/r)^(2*(j-1)+3))*(2*(j-1)+1);
        end
    end
elseif var==2
    if r<=a
        for j=2:n
            Q=P;
            P=Q + coeff(2*j-1)*Pol(2*j-1)*((r/a)^(2*(j-1)-2))*2*(j-1)/(2*(j-1)-1);
        end
    else
        for j=2:n
            Q=P;
            P=Q + coeff(2*j-1)*Pol(2*j-1)*((a/r)^(2*(j-1)+3))*(2*(j-1)+1)/j;
        end
    end
else
    if r<=a
        for j=2:n
            Q=P;
            P=Q + coeff(2*j-1)*Pol(2*j-1)*((r/a)^(2*(j-1)-2))*2*(j-1)/(2*(j-1)-1);
        end
    else
        for j=1:n
            Q=P;
            P=Q + coeff(2*j-1)*Pol(2*j-1)*((a/r)^(2*(j-1)+3))*(2*(j-1)+1)/(2*(j-1)+2);
        end
    end
end
end

function dP = Pol_Leg_Derivative_Sum(r,a,n,var,theta,coeff,dPol)
dP=0;
if var==3
    if r<=a
        for j=2:n
            Q=dP;
            dP=Q + coeff(2*j-1)*(cos(theta)/sin(theta))*dPol(2*j-1)*((r/a)^(2*(j-1)-2))*(1/(2*(j-1)-1));
        end
    else
        for j=1:n
            Q=dP;
            dP=Q + coeff(2*j-1)*(cos(theta)/sin(theta))*dPol(2*j-1)*((a/r)^(2*(j-1)+3))*(1/(2*(j-1)+2));
        end
    end
elseif var==4
    if r<=a
        for j=2:n
            Q=dP;
            dP=Q + coeff(2*j-1)*dPol(2*j-1)*((r/a)^(2*(j-1)-2));
        end
    else
        for j=2:n
            Q=dP;
            dP=Q + coeff(2*j-1)*dPol(2*j-1)*((a/r)^(2*(j-1)+3));
        end
    end
end
end

function ddP = Pol_Leg_2Derivative_Sum(r,a,n,var,coeff,ddPol)
ddP=0;
if var==2
    if r<=a
        for j=2:n
            Q=ddP;
            ddP=Q + coeff(2*j-1)*ddPol(2*j-1)*((r/a)^(2*(j-1)-2))*(1/(2*(j-1)-1));
        end
    else
        for j=2:n
            Q=ddP;
            ddP=Q + coeff(2*j-1)*ddPol(2*j-1)*((a/r)^(2*(j-1)+3))*(1/j);
        end
    end
else
    ddP=0;
end
end
