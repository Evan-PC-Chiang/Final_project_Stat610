
%elastic thickness
H= 20;

%upper ramp
xtip = 0;
ztip = 0;
dip = 30;
zbot = 5;


%middle flat
dip_middle = 5;
zbot_middle = 10;


%lower ramp
dip_lower=30;

xobs = linspace(-300,300,200);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




mu=1; t=1; zobs=0:H; 



%upper ramp
slip = 1;
L = (zbot-ztip)/sin(dip*pi/180);
[ux1,uz1,ux1e,uz1e,x1] = drive_visc_ramp_body_return_elastic(H,ztip,0,dip,L,slip,mu,t,zobs);
x1 = x1+xtip;

upper_ends = [xtip ztip; xtip+L*cos(dip*pi/180) ztip+L*sin(dip*pi/180) ];

%middle ramp
slip = 1;
xtip = xtip + L*cos(dip*pi/180);
L = (zbot_middle-zbot)/sin(dip_middle*pi/180);
ztip = zbot;
[ux2,uz2,ux2e,uz2e,x2] = drive_visc_ramp_body_return_elastic(H,ztip,0,dip_middle,L,slip,mu,t,zobs);
x2 = x2+xtip;

middle_ends = [xtip ztip; xtip+L*cos(dip_middle*pi/180) ztip+L*sin(dip_middle*pi/180) ];


%lower ramp
slip = 1;
xtip = xtip + L*cos(dip_middle*pi/180);
L = (H-zbot_middle)/sin(dip_lower*pi/180);
ztip = zbot_middle;
[ux3,uz3,ux3e,uz3e,x3] = drive_visc_ramp_body_return_elastic(H,ztip,0,dip_lower,L,slip,mu,t,zobs);
x3 = x3+xtip;

lower_ends = [xtip ztip; xtip+L*cos(dip_lower*pi/180) ztip+L*sin(dip_lower*pi/180) ];


%interplate to xobs and add
[Xobs,Zobs] = meshgrid(xobs,zobs);
[X1,Z1] = meshgrid(x1,zobs);
[X2,Z2] = meshgrid(x2,zobs);
[X3,Z3] = meshgrid(x3,zobs);
ux1 = griddata(X1,Z1,ux1,Xobs,Zobs);%,'nearest');
uz1 = griddata(X1,Z1,uz1,Xobs,Zobs);%,'nearest');
ux2 = griddata(X2,Z2,ux2,Xobs,Zobs);%,'nearest');
uz2 = griddata(X2,Z2,uz2,Xobs,Zobs);%,'nearest');
ux3 = griddata(X3,Z3,ux3,Xobs,Zobs);%,'nearest');
uz3 = griddata(X3,Z3,uz3,Xobs,Zobs);%,'nearest');

ux1e = griddata(X1,Z1,ux1e,Xobs,Zobs);%,'nearest');
uz1e = griddata(X1,Z1,uz1e,Xobs,Zobs);%,'nearest');
ux2e = griddata(X2,Z2,ux2e,Xobs,Zobs);%,'nearest');
uz2e = griddata(X2,Z2,uz2e,Xobs,Zobs);%,'nearest');
ux3e = griddata(X3,Z3,ux3e,Xobs,Zobs);%,'nearest');
uz3e = griddata(X3,Z3,uz3e,Xobs,Zobs);%,'nearest');

ux = ux1 + ux2 + ux3;
uz = uz1 + uz2 + uz3;

uxe = ux1e + ux2e + ux3e;
uze = uz1e + uz2e + uz3e;

figure
quiver(Xobs,-Zobs,ux,uz)
hold on
plot(upper_ends(:,1),-upper_ends(:,2),'r')
plot(middle_ends(:,1),-middle_ends(:,2),'r')
plot(lower_ends(:,1),-lower_ends(:,2),'r')
axis equal


figure
subplot(211)
plot(xobs,ux(1,:))
grid on
subplot(212)
plot(xobs,uz(1,:))
grid on



figure
quiver(Xobs,-Zobs,uxe,uze)
hold on
plot(upper_ends(:,1),-upper_ends(:,2),'r')
plot(middle_ends(:,1),-middle_ends(:,2),'r')
plot(lower_ends(:,1),-lower_ends(:,2),'r')
axis equal


figure
subplot(211)
plot(xobs,uxe(1,:))
grid on
subplot(212)
plot(xobs,uze(1,:))
grid on
