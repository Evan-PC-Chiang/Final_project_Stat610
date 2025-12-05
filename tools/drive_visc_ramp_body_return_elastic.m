function [u1,u2,edgu1,edgu2,x] = drive_visc_ramp_body_return_elastic(H,ztip,xtip,dip,L,slip,mu,t,zobs)

tR=1;
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
   
  
[U1(n,:,:),U2(n,:,:), x] = visc_dip_ramp_body(H,xs(n),zs(n),M,t,tR,N,mu,slip,dip,zobs);

end %n (loop over number of point sources)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%integrate along fault
wf=repmat(wf,[1 1 length(zobs)]);
u1=sum(wf.*U1,1);      %wf -- weighting factors obtained above
u2=sum(wf.*U2,1);
u1=squeeze(u1)';u2=squeeze(u2)';


%calculate elastic part of solution

x1=x;
x2=-zobs;
[x1,x2]=meshgrid(x1,x2);

m=[xtip -ztip dip*180/pi slip*t];
[edgu1_1, edgu2_1] =EdgeDisp(m,x1,x2,0.25);

m=[xtip+L*cos(dip) -(ztip+L*sin(dip)) dip*180/pi -slip*t];
[edgu1_2, edgu2_2] =EdgeDisp(m,x1,x2,0.25);

edgu1=edgu1_1+edgu1_2;
edgu2=edgu2_1+edgu2_2;

u1=u1+edgu1;
u2=u2+edgu2;

