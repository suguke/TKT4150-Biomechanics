clear
%--------------------- VARIABLES/CONSTANTS -------------------------%
a=0.5e-2;   mu=1.05e-3; T=0.5;  f=1/T;  w=2*pi*f;   rho=1e3; 
Ns = 100;                       %-Radial resolution
Nt = 200;                        %-Number of time levels
Np = 4;                        
Tmin=0.0;   Tmax=Np*T;
                         
r = linspace(0,a,Ns);           %-Make radial vector r
t = linspace(Tmin,Tmax,Nt);     %-Make time vector t

alpha = sqrt(rho*w/mu)*a      %-Womersley number
p0=1.0;
p=p0*sin(w*t).^2;               %-Make a time varying complex pressure 
                                % vector dp/dz(w,t), use a constant 
                                % amplitude p0=1.0
dp = p0*exp(1i*w*t);
                                
%------ CALCULATE WOMERSLEY PROFILES WITH BESSEL FUNCTIONS ----------%
for i=1:Nt
%v(i,:)=...%                    %-Use the built in matlab function
                                % besselj(nu,Z) to calculate the Bessel
                                % functions needed to find an expression 
                                % for the velocity vector v(i,:) as a 
                                % function of the pressure function p(i). 
    alpha=2
    v2(i,:) = (1i/w/rho).*dp(i).*(1-besselj(0,1i^(1.5).*alpha.*r./a)./besselj(0,1i^(1.5).*alpha));
    alpha=5
    v5(i,:) = (1i/w/rho).*dp(i).*(1-besselj(0,1i^(1.5).*alpha.*r./a)./besselj(0,1i^(1.5).*alpha));
    alpha=10
    v10(i,:) = (1i/w/rho).*dp(i).*(1-besselj(0,1i^(1.5).*alpha.*r./a)./besselj(0,1i^(1.5).*alpha));
    alpha=20
    v20(i,:) = (1i/w/rho).*dp(i).*(1-besselj(0,1i^(1.5).*alpha.*r./a)./besselj(0,1i^(1.5).*alpha));
end


%-------------------- PLOT WOMERSLEY PROFILES -----------------------%
p    = -real(p)/p0;
vmax = max(max(real(v2)));
v2    = real(v2)/vmax;

vmax = max(max(real(v5)));
v5    = real(v5)/vmax;
vmax = max(max(real(v10)));
v10    = real(v10)/vmax;
vmax = max(max(real(v20)));
v20    = real(v20)/vmax;
r = [-r(Ns:-1:1) r];            %-Make r a vector from -a to a
v2 = [v2(:,Ns:-1:1) v2(:,:)];      %-Mirror v about r=0
v5 = [v5(:,Ns:-1:1) v5(:,:)];  
v10 = [v10(:,Ns:-1:1) v10(:,:)];  
v20 = [v20(:,Ns:-1:1) v20(:,:)];  

onesArray = ones(size(r));
timeLevels = 1:Nt;
% figure()
% for timeLevel=timeLevels
%     plot(r/a,v(timeLevel,:)+timeLevel-1,'k');
%     plot(r/a,timeLevel*ones(size(r))-1, 'k')
%     hold all
%     idx = 1:10:2*Ns;
%     quiver(r(idx)/a,timeLevel-1,0*r(idx),v(timeLevel,idx),1)
%     hold all
% end

frames = moviein(length(timeLevels));
for i=1:length(timeLevels)
    timeLevel = timeLevels(i)
    clf
    plot(r/a,v2(timeLevel,:))
    hold all
    plot(r/a,v5(timeLevel,:))
    plot(r/a,v10(timeLevel,:))
    plot(r/a,v20(timeLevel,:))
    set(gca,'ylim',[-1 1]);
    legend('\alpha=2','\alpha=5','\alpha=10','\alpha=20')
%     idx = 1:10:2*Ns;
%     quiver(r(idx)/a,0,0*r(idx),v(timeLevel,idx),1);set(gca,'ylim',[-1 1]);
%    tstr = sprintf('\\alpha=%f, t=%f',alpha, t(timeLevel))
    tstr = sprintf('t=%f', t(timeLevel))
    title(tstr)
    frames(:,i) = getframe; 
end
% %%
% figure
% xlabel('r')
% ylabel('v')
% movie(frames,T*Np,Np*T/Nt)



