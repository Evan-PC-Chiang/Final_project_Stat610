function [u1,u2,x] = drive_visc_cycle_body(H,ztip,xtip,dip,L,slip,mu,t,T,zobs)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[u1,u2] = drive_visc_cycle(H,ztip,xtip,dip,L,slip,mu,t,tR,T)
%
%%%input parameters:
%H             depth (POSITIVE) to bottom of layer1 (km)
%ztip          depth (POSITIVE) of fault tip (km)
%xtip          horizontal position of fault tip (km)
%dip           dip of fault (degrees)
%L             length of fault (km)
%slip          amount of slip (m)
%mu				shear modulus
%t  		      time normalized by tR (tR = relaxation time)
%T             recurrance interval normalized by tR
%
%%%Outputs:
%u1,u2         horizontal and vertical displacements at 500 points
%              spaced evenly from -20H to 20H
% Kaj Johnson, June 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
[U1(n,:,:),U2(n,:,:), x] = visc_dip_cycle_body(H,xs(n),zs(n),M,T,t,tR,N,mu,zobs);
end %n (loop over number of point sources)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%integrate along fault
wf=repmat(wf,[1 1 length(zobs)]);
u1=sum(wf.*U1,1);      %wf -- weighting factors obtained above
u2=sum(wf.*U2,1);
u1=squeeze(u1)';u2=squeeze(u2)';

