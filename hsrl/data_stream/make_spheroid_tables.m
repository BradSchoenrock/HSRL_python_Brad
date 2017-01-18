function [spheroid_tables]=make_spheroid_tables(temperature,rho_air,spheroid_prams,fig_offset,phase);
%[spheroid_tables]=make_spheroid_tables(temperature,rho_air,spheroid_prams,fig_offset,phase);
%size distribution N(D)dD=N*k*D^alpha*exp(-bD^gamma)dD
%                       k=gamma*b^((alpha+1)/gamma)/Gamma((alpha+1)/gamma)
%spheroids height described by power law h=sigma_h*D_r*(D/D_r)^eta
%spheroid_prams          = structure describing particles
%phase                   = 'ice' or 'water'                  
%spheroid_prams.zeta     = power 
%              .sigma_h  = ratio of spheroid height to diameter ad D=D_r
%              .Dr       = reference diameter used to make sigma_h non-dimensional
%              .alpha_ice= parameter in Gamma distribution for ice particle sizes
%              .alpha_h20= parameter in Gamma distribution for water droplets
%              .gamma_ice= parameter in Gamma distribution for ice particles
%              .gamma_h20= parameter in Gamma distribution for water dropletes
%spheroid_tables.deff_prime             = deff_prime vector (m), log spacing              
%		.mode_diameter          = mode_diameter vector (m), log spacing
%		.zeta                   = power law vector,  Auer and Veal exponent
%		.Vrf_vs_deff_prime_zeta (m/s), size(length(length(zeta),length(deff_prime))
%		.zeta_vs_vrf_deff_prime (non-dim), size(length(Vrf),length(deff_prime))
%		.spheroid_prams



%constants, program works in cgs with input and output in mks
rho_ice=0.91; %gr/cm^3
rho_air=rho_air*1e-3; %convert kg/m^3 to gr/cm^3 

if strcmp(phase,'ice')
  delta_0=5.83;
  C_0=0.6; %page 4345 Khvorostyanov and Curry, Dec 2005 JAS
  rho=rho_ice;
  gam=spheroid_prams.gamma_ice;
  alpha=spheroid_prams.alpha_ice;
elseif strcmp(phase,'water')
  rho=1;
  delta_0=9.06;
  C_0=0.292;
  gam=spheroid_prams.gamma_h20;
  alpha=spheroid_prams.alpha_h20;
end

C_1=4/(sqrt(C_0)*delta_0^2);
C_2=0.25*delta_0^2;
g=980;  %acceleration of gravity in cgs
a_0=1.7e-3; %page 4345 K&C
b_0=0.8;

%dynamic viscosity in poise, gr/(cm sec), K&C page 4347
if temperature-273 < 0
  eta=1.718e-4*(1+0.00285*(temperature-273)-6.9e-6*(temperature-273)^2);
else
  eta=1.718e-4*(1+0.00285*(temperature-273));
end

spheroid_prams
%particle parameters
sigma_h  =spheroid_prams.sigma_h;
Dr       =spheroid_prams.Dr*1e2; %convert meters to cm 


%effective diameter prime vector
deff_prime=logspace(0,4,100)*1e-4; % 1 micron to 1 cm expressed in cm



n_zeta=10;
for i=1:n_zeta+1  %loop over zeta values between 0 and 1
  zeta(i)=(i-1)/n_zeta;

  gamma_ratio=gamma((alpha+3)/gam)/gamma((2*zeta(i)+alpha+5)/gam);

  %mode diameters corresponding to the above deff_primes
  Dm=sigma_h^(-1/(zeta(i)+1))*Dr^((zeta(i)-1)/(zeta(i)+1))*(alpha/gam)^(1/gam)...
     *gamma_ratio^(1/(2*zeta(i)+2))*deff_prime.^(2/(zeta(i)+1));

  if 0
    figure(1)
    loglog(deff_prime,Dm)
    title('Mode diameter vs. deff prime')
    xlabel('Deff prime (cm)')
    grid on
  end

  %D is a vector of diameters used for the numeric integration
  D=logspace(-1,4,500)*1e-4; %D from 0.1 micron to 1 cm

  % X=best number vector
  X_const=4*sigma_h*Dr^(1-zeta(i))*rho*rho_air*g/(3*eta^2);
  X=X_const*D.^(2+zeta(i));

  if 0
    figure(2)
    loglog(D,X)
    title('Best number vs Diameter')
    xlabel('Diameter (cm)')
    grid on
  end

  %fall velocities for single particles
  a_0=0;
  Vf(:,i)=(eta./(rho_air*D)).*(C_2*(((1+C_1.*X.^0.5).^0.5)-1).^2-a_0*X.^b_0);

  %compute radar weighted and mass weighted fall velocity
  for k=1:length(Dm)
    %radar weighted fall velocity
    Vrf_table(i,k)=sum(Vf(:,i)'.*exp(-(alpha/gam).*(D/Dm(k)).^gam).* ...
	     D.^(4+2*zeta(i)+alpha));
    Vrf_table(i,k)=Vrf_table(i,k)./(sum(exp(-(alpha/gam).*(D/Dm(k)).^gam).* ...
	     D.^(4+2*zeta(i)+alpha)));
    %mass weighted fall velocity
    Vmf_table(i,k)=sum(Vf(:,i)'.*exp(-(alpha/gam).*(D/Dm(k)).^gam).* ...
	     D.^(2+zeta(i)+alpha));
    Vmf_table(i,k)=Vrf_table(i,k)./(sum(exp(-(alpha/gam).*(D/Dm(k)).^gam).* ...
	     D.^(2+zeta(i)+alpha)));
  end
end

%make table of zeta vs radar_fall_velocity and deff_prime
%and  table of mass_weighted_fall_velocity/radar_weighted_fall_velocity
%    vs vrf and deff_prime

Vrf=0:150;
Vrf=Vrf*2;  % radar fall velocity 2 cm/s to 300 cm/sec

%Vrf_table(i+1,:)=Vrf(length(Vrf));
Vrf_table(1,:)=0;
for k=1:length(deff_prime)
  zeta_vs_vrf_deff_prime(:,k)=interp1(Vrf_table(:,k),zeta,Vrf,'linear',1);
  Vmf_Vrf_vs_vrf_deff_prime(:,k)=interp1(Vmf_table(:,k)./Vrf_table(:,k),zeta,Vrf,'linear',1);
end

j=1;
legend_str=['legend(''',num2str(zeta(1))];
while j<length(zeta)
  j=j+1;
  legend_str=[legend_str,'''',',','''', num2str(zeta(j))];
end
legend_str=[legend_str,'''',',''location'',''northwest'')'];


if 0
 figure(300)
 loglog(D*1e4,Vf/100)
 axis([10 1e4 .01 2])
 title('fall velocity vs Diameter')
 xlabel('Diameter (microns)')
 ylabel('Fall velocity (m/s)')
 grid on
 
 figure(400)
 semilogx(D*1e4,Vf/100)
 axis([10 1e4 0 20])
 grid on
 title('fall velocity vs Diameter')
 xlabel('Diameter (micron)')
 ylabel('Fall velocity (m/s)')
 eval(legend_str)
end

if ~isempty(fig_offset)
figure(500+fig_offset)
semilogx(deff_prime*1e4,Vrf_table/100)
axis([1 1e3 -0.5 3])
grid
xlabel('D_e_f_f'' (microns)')
ylabel('Radar weighted fall velocity (m/s)')
      
eval(legend_str);
end

if 0
figure(600)
semilogx(Dm*1e4,Vrf_table/100)
axis([10 1e4 0 3])
grid
xlabel('Mode diameter (microns)')
ylabel('Radar weighted fall velocity (m/s)')
end      




%convert return values to mks
spheroid_tables=struct('deff_prime',deff_prime/100 ...
		       ,'mode_diameter',Dm/100 ...
		       ,'Vrf',Vrf/100 ...
		       ,'zeta',zeta...
		       ,'Vrf_vs_deff_prime_and_zeta',Vrf_table/100 ...
		       ,'spheroid_prams',spheroid_prams ...
                       ,'zeta_vs_vrf_and_deff_prime',zeta_vs_vrf_deff_prime ...
		       ,'Vmf_Vrf_vs_vrf_and_deff_prime',Vmf_Vrf_vs_vrf_deff_prime);
