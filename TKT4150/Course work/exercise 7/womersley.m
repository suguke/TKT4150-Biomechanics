%--------------------- VARIABLES/CONSTANTS -------------------------%
a=0.5e-2;   mu=1.05e-3; T=0.5;  f=1/T;  w=2*pi*f;   rho=1e3; 
Ns = 100;                       %-Radial resolution
Nt = 50;                        %-Number of time levels
Np = 10;                        
Tmin=0.0;   Tmax=Np*T;
                         
r = linspace(0,a,Ns);           %-Make radial vector r
t = linspace(Tmin,Tmax,Nt);     %-Make time vector t

alpha = sqrt(rho*w/mu)*a;       %-Womersley number
p0=1.0;
p=p0*sin(w*t).^2%               %-Make a time varying complex pressure 
                                % vector dp/dz(w,t), use a constant 
                                % amplitude p0=1.0

                                
%------ CALCULATE WOMERSLEY PROFILES WITH BESSEL FUNCTIONS ----------%
for i=1:Nt
%v(i,:)=...%                    %-Use the built in matlab function
                                % besselj(nu,Z) to calculate the Bessel
                                % functions needed to find an expression 
                                % for the velocity vector v(i,:) as a 
                                % function of the pressure function p(i). 
end


%-------------------- PLOT WOMERSLEY PROFILES -----------------------%
p    = -real(p)/p0;
vmax = max(max(real(v)));
v    = real(v)/vmax;

r = [-r(Ns:-1:1) r];            %-Make r a vector from -a to a
v = [v(:,Ns:-1:1) v(:,:)];      %-Mirror v about r=0

onesArray = ones(size(r));
timeLevel = 1;

subplot(1,2,1);h(1)=plot(r/a,onesArray*p(timeLevel),'r');set(gca,'ylim',[-1 1]);
subplot(1,2,2);h(2)=plot(r/a,v(timeLevel,:));set(gca,'ylim',[-1 1]);
set(h(:),'linewidth',3);




