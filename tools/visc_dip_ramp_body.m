	function [U1,U2, x] = visc_dip_ramp_body(H,xs,zs,MT,t,tR,N,mu,slip,dip,zobs);
   
   
% [U1,U2, x] = visc_dip(H,xs,zs,M,t,tR,N);
%
% INPUT:
%	H	= Layer thickness
%	xs	= x position of point source
%	zs	= Depth of point source, < H.
%	M	= moment tensor.
%	t   = time since earthquake
%	tR  = relaxation time 
%	N   = Number of points in FFT
%
% OUTPUT:
%	u1	= displacement
%	u2	= displacement
%	x	= x-coordinate corresponding to u
%
% Propagator matrix method for viscoelastic 2-D dip-slip faults
% This method uses real valued A
% Paul Segall, sept. 2001

xmax = 20*H;       %range of x
delta = 2*xmax/N;  %sampling interval in x
kmax = 0.5/delta;  %maximum wavenumber
x = linspace(-xmax, xmax,N);
ks = linspace(-kmax,kmax,N)*(2*pi);  
% factor of 2*pi accounts for different definition of transform

% some material constants
nu = 0.25;
lam=-2*nu.*mu./(2*nu-1);
K=lam+2*mu/3;   %bulk modulus
g=lam+2*mu;
B = 1/tR;

u1 = zeros(size(ks));
u2 = zeros(size(ks));

rg=3*10^-3;
G=[1 0 0 0;0 1 0 0;0 0 1 0;0 -rg 0 1];

for j=1:N
   
	k = ks(j);
   
	%set up A matrix in equation dv=A*v+f
         A=[0  k 1/mu 0; ...
        -k*lam/g 0 0 1/g; ...
        4*k^2*mu*(lam+mu)/g 0 0 k*lam/g; ...
        0 0 -k 0];
  
  
  	%propagator from top of halfspace to surface
	z=0; z0=H;
	Ph = G*expm(A*(z-z0));
   
   %propagator from source to surface
	z=0; z0=zs;
	Pzs= G*expm(A*(z-z0));
   
   
	%body forces
	F = [i*MT(1,2)/mu; MT(2,2)/g; i*i*k*MT(1,1) + k*lam*MT(2,2)/g; 0]*exp(-i*xs*k);
   
   PF=Pzs*F;
   
	%compute roots
	[s1, s2, s3, a] = roots_right(k,mu,K,lam,B,Ph);
	%[s1, s2, s3, a] = roots_kaj(k,mu,K,B,Ph);
   
   t0=1000; t1=1000+t;  %start sliding 1000 relaxation times ago -- reach steady state
   
	%% RETURN VISCOELASTIC PART ONLY
	%% FIRST root:
	V = makeV_right(k,K,lam,s1,B,H,mu);
   M = Ph(3:4,:)*V;
	Adj(1,1) = M(2,2); Adj(2,2) = M(1,1); 
		Adj(1,2) = -M(1,2); Adj(2,1) = -M(2,1);
	Num = -Ph*V*Adj*Pzs(3:4,:)*F;
	Num1 = -V*Adj*Pzs(3:4,:)*F;
	Den1 = a*s1^2*(s1-s2)*(s1-s3);
  term1 = Num*(exp(s1*t1)-exp(s1*t0))/Den1;
  
	%% SECOND root:
	V = makeV_right(k,K,lam,s2,B,H,mu);
	M = Ph(3:4,:)*V;
	Adj(1,1) = M(2,2); Adj(2,2) = M(1,1); 
		Adj(1,2) = -M(1,2); Adj(2,1) = -M(2,1);
	Num = -Ph*V*Adj*Pzs(3:4,:)*F;
	Num2 = -V*Adj*Pzs(3:4,:)*F;
	Den2 = a*s2^2*(s2-s1)*(s2-s3);
	term2 = Num*(exp(s2*t1)-exp(s2*t0))/Den2;
   
   
   %% THIRD root:
   V = makeV_right(k,K,lam,s3,B,H,mu);
	M = Ph(3:4,:)*V;
	Adj(1,1) = M(2,2); Adj(2,2) = M(1,1); 
		Adj(1,2) = -M(1,2); Adj(2,1) = -M(2,1);
	Num = -Ph*V*Adj*Pzs(3:4,:)*F;
	Num3 = -V*Adj*Pzs(3:4,:)*F;
	Den3 =  a*s3^2*(s3-s1)*(s3-s2);
	term3 = Num*(exp(s3*t1)-exp(s3*t0))/Den3;
   
       %%zero root
   V = makeV_right(k,K,lam,0,B,H,mu);
	M = Ph(3:4,:)*V;
	Adj(1,1) = M(2,2); Adj(2,2) = M(1,1); 
		Adj(1,2) = -M(1,2); Adj(2,1) = -M(2,1);
	Num = -Ph*V*Adj*Pzs(3:4,:)*F;
	Num5 = -V*Adj*Pzs(3:4,:)*F;
	Den5 =  a*(-s1)*(-s2)*(-s3);
	term5 = Num/Den5;

            %%zero root (order 2 zero root)
			if k<0
   			term4 = res1_lt0(k,K,lam,0,B,H,mu,PF,Ph,s1,s2,s3,a);
			else
   			term4 = res1_gt0(k,K,lam,0,B,H,mu,PF,Ph,s1,s2,s3,a);
			end

   
Fhat = term1 + term2 + term3 + PF*t + term5*t;



%%%elastic solution
 A=[0 k 1/mu(1) 0;-k*lam(1)/g(1) 0 0 1/g(1);4*k^2*mu(1)*(lam(1)+mu(1))/g(1) 0 0 k*lam(1)/g(1);0 0 -k 0];
%propagator from top of halfspace to surface
z=0;
z0=H;
Ph = G*expm(A*(z-z0));

%propagator from source to surface
z=0;
z0=zs;
Pzs = G*expm(A*(z-z0));


%basis vectors for solution to homogeneous equation
if k<0
   a1=[-.5/mu .5/mu -k k]';
   b1=[-(2*lam+3*mu)/(2*mu*k*(lam+mu));lam/(2*mu*k*(lam+mu));-2;1];
   d1=a1;
   d2=b1+a1*H;
else
   a2=[-.5/mu -.5/mu k k]';
   b2=[1/(2*k*(lam+mu));-(lam+2*mu)/(2*mu*k*(lam+mu));0;1];
   d1=a2;
   d2=b2+a2*H;
end

       
%body forces
M11=(slip*cos(dip)*sin(dip)*2*mu)*t;
M12=(slip*sin(dip)*sin(dip)*mu-slip*cos(dip)*cos(dip)*mu)*t;
M22=(-slip*sin(dip)*cos(dip)*2*mu)*t;

f=[i*M12/mu;M22/g;-k*M11+k*M22*lam/g;0]*exp(-i*xs*k);

Pf=Pzs*f;
%calculate constants using traction free boundary condition
Pd1=Ph*d1;
Pd2=Ph*d2;
M=[Pd1(3) Pd2(3);Pd1(4) Pd2(4)];
b=-Pf(3:4);
unknown=inv(M)*b;
c1=unknown(1); %c1,c2 contstants
c2=unknown(2);

PV=c1*Pd1+c2*Pd2;
Fhate=PV+Pf;
Fobse=c1*d1+c2*d2;

%propagate solution to zobs
   
   	%set up A matrix in equation dv=A*v+f
         A=[0  k 1/mu 0; ...
        -k*lam/g 0 0 1/g; ...
        4*k^2*mu*(lam+mu)/g 0 0 k*lam/g; ...
        0 0 -k 0];
  
for loop=1:length(zobs);
   if zobs(loop)<zs
      	z=zobs(loop); z0=0;
         Pzobs = expm(A*(z-z0))*inv(G);
			Fobs=Pzobs*(Fhat - Fhate); %remove elastic part
         Fstress=Pzobs*Fhat;
         u1(loop,j) = -i*Fobs(1);
         u2(loop,j) = -Fobs(2);
      else
         z=zobs(loop); z0=H;
         Pzobs = expm(A*(z-z0));
         term1 = Pzobs*Num1*(exp(s1*t1)-exp(s1*t0))/Den1;
         term2 = Pzobs*Num2*(exp(s2*t1)-exp(s2*t0))/Den2;
         term3 = Pzobs*Num3*(exp(s3*t1)-exp(s3*t0))/Den3;
       	term5 = Pzobs*Num5/Den5;
         Fobs = term1 + term2 + term3 + term5*t - Pzobs*Fobse; %remove elastic part
         u1(loop,j) = -i*Fobs(1);
         u2(loop,j) = -Fobs(2);
      end
end %loop


end %j (loop over wave number)


for loop=1:length(zobs)

%filter out high frequency noise at large k
temp=real(u1(loop,:));
tempo=mean(abs(temp(end-10:end)));
temp1=zeros(size(temp));temp2=zeros(size(temp));temp3=zeros(size(temp));temp4=zeros(size(temp));temp5=zeros(size(temp));
if tempo>1e-3
   temp=diff(temp,10);
   temp1=temp(1:5:end);temp2=temp(2:5:end);temp3=temp(3:5:end);temp4=temp(4:5:end);temp5=temp(5:5:end);
   temp=temp1(1:N/10)+temp2(1:N/10)+temp3(1:N/10)+temp4(1:N/10)+temp5(1:N/10);
   
   index=abs(temp)<10^-4;
   
  
   
   q=0;
   p=0;
   while q==0
      p=p+1;
      q=index(p);
   end

u1(loop,1:5*p)=0;
u1(loop,end-5*p:end)=0;
u2(loop,1:5*p)=0;
u2(loop,end-5*p:end)=0;

end



      %inverse FT into physical space

    %divide by sampling interval to turn analytical transform
    %shift negative frequencies to N/2+1 to N-1
    
    %inverse FFT
	u1s=ifft(fftshift(u1(loop,:)/delta));
	u2s=ifft(fftshift(u2(loop,:)/delta));
   
    %shift the solution back
	U1(1,:,loop)=real(fftshift(u1s));
	U2(1,:,loop)=real(fftshift(u2s));
end

  


