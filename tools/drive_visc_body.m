%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input parameters:
H=1;             %depth (POSITIVE) to bottom of layer1 (km)
ztip=0.5;          %depth (POSITIVE) of fault tip (km)
xtip=0;          %horizontal position of fault tip (km)
dip=30;          %dip of fault (degrees)
L=0.5*H/sin(dip*pi/180);             %length of fault (km)
slip=1;          %amount of slip (m)
mu = 1;

t = 2;  %time normalized by tR
tR = 1;		% relaxation time

dip=dip*pi/180;
nps=10;          %number of point sources
D=ztip+L*sin(dip);

N = 1000;


%setup spacing of points on fault for Gaussian integration
u=(1:nps-1)./sqrt((2*(1:nps-1)).^2-1);
[vc,bp]=eig(diag(u,-1)+diag(u,1));
[bp,K]=sort(diag(bp));
a=0;b=L;
wf=(2*vc(1,K).^2)*(b-a)/2;      %weighting factors used in the summation
wf=wf';
%wf=repmat(wf,1,size(U1,2));
wf=repmat(wf,1,N);
xp=(a+b)/2+(b-a)/2*bp;  %spacing of integration points (point sources)


%% set up moment tensor
M = zeros(2,2);
mus = mu;
M(1,1)=slip*cos(dip)*sin(dip)*2*mus;
M(1,2)=slip*sin(dip)*sin(dip)*mus-slip*cos(dip)*cos(dip)*mus;
M(2,2)=-slip*sin(dip)*cos(dip)*2*mus;


%for kk = 1:length(t)   %%loop on time



%loop of number of point sources (nps)
for n=1:nps
   xs(n)=xtip+xp(n)*cos(dip);
   zs(n)=ztip+xp(n)*sin(dip);              %position of point source (xs,zs)

%[U1(n,:),U2(n,:), x] = visc_dip_grav2(H,xs(n),zs(n),M,t,tR,N,mu);
zobs=0.1:0.1:1;
[U1(n,:,:),U2(n,:,:), x] = visc_dip_stresses(H,xs(n),zs(n),M,t,tR,N,mu,zobs);
%[U2(n,:), x] = visc_dip_numerical(slip,dip,H,xs(n),zs(n),M,t,tR,N,mu);
%[U1(n,:),U2(n,:), x] = visc_dip_ramp(H,xs(n),zs(n),M,t,tR,N,mu);
%[U1(n,:), x] = visc_dip_ramp_surf(slip,dip,H,xs(n),zs(n),M,t,tR,N,mu,L);

end %n (loop over number of point sources)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%integrate along fault

%UU1=sum(wf.*U1,1);      %wf -- weighting factors obtained above
%UU2=sum(wf.*U2,1);

wf=repmat(wf,[1 1 length(zobs)]);


u1=sum(wf.*U1,1);      %wf -- weighting factors obtained above
u2=sum(wf.*U2,1);

u1=squeeze(u1)';u2=squeeze(u2)';
