
%elastic thickness
H= 20;

%upper ramp
xtip_up = 0;
ztip_up = 0;
dip_up = 30;
zbot_up = 5;


%middle flat
dip_middle = 5;
zbot_middle = 10;


%lower ramp
dip_lower=30;

%slip rate
srate = 1;

xobs = linspace(-300,300,200);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




mu=1; t=1; zobs=0:H; 



%% upper ramp
slip = srate;
L_up = (zbot_up-ztip_up)/sin(dip_up*pi/180);
[ux1,uz1,ux1e,uz1e,x1] = drive_visc_ramp_body_return_elastic(H,ztip_up,0,dip_up,L_up,slip,mu,t,zobs); %shift tip of fault to x=0
x1 = x1+xtip_up; 

upper_ends = [xtip_up ztip_up; xtip_up+L_up*cos(dip_up*pi/180) ztip_up+L_up*sin(dip_up*pi/180) ];

%% middle ramp
slip = srate;
xtip_middle = xtip_up + L_up*cos(dip_up*pi/180);
L_middle = (zbot_middle-zbot_up)/sin(dip_middle*pi/180);
ztip_middle = zbot_up;
[ux2,uz2,ux2e,uz2e,x2] = drive_visc_ramp_body_return_elastic(H,ztip_middle,0,dip_middle,L_middle,slip,mu,t,zobs);  %shift tip of fault to x=0
x2 = x2+xtip_middle;

middle_ends = [xtip_middle ztip_middle; xtip_middle+L_middle*cos(dip_middle*pi/180) ztip_middle+L_middle*sin(dip_middle*pi/180) ];


%% lower ramp
slip = srate;
xtip_lower = xtip_middle + L_middle*cos(dip_middle*pi/180);
L_lower = (H-zbot_middle)/sin(dip_lower*pi/180);
ztip_lower = zbot_middle;
[ux3,uz3,ux3e,uz3e,x3] = drive_visc_ramp_body_return_elastic(H,ztip_lower,0,dip_lower,L_lower,slip,mu,t,zobs); %shift tip of fault to x=0
x3 = x3+xtip_lower;

lower_ends = [xtip_lower ztip_lower; xtip_lower+L_lower*cos(dip_lower*pi/180) ztip_lower+L_lower*sin(dip_lower*pi/180) ];


%% axial shear zones

%upper axial shear zone
dip_axial_up =  90 - (dip_middle + dip_up)/2;
slip = 2*srate*sin(abs(dip_up-dip_middle)/2*pi/180);
xtip_axial_up = xtip_middle + ztip_middle/tan(dip_axial_up*pi/180); 
L_axial_up = ztip_middle/sin(dip_axial_up*pi/180);
[ux4,uz4,ux4e,uz4e,x4] = drive_visc_ramp_body_return_elastic(H,0,0,180-dip_axial_up,L_axial_up,-slip,mu,t,zobs); %shift tip of fault to x=0, dip-180 dip and negative sign on slip for surface dipping down to left
x4 = x4+xtip_axial_up;

axial_up_ends = [xtip_axial_up 0; xtip_middle ztip_middle ];


%lower axial shear zone
dip_axial_lower =  90 - (dip_middle + dip_lower)/2;
slip = -2*srate*sin(abs(dip_middle-dip_lower)/2*pi/180);
xtip_axial_lower = xtip_lower + ztip_lower/tan(dip_axial_lower*pi/180); 
L_axial_lower = ztip_lower/sin(dip_axial_lower*pi/180);
[ux5,uz5,ux5e,uz5e,x5] = drive_visc_ramp_body_return_elastic(H,0,0,180-dip_axial_lower,L_axial_lower,-slip,mu,t,zobs); %shift tip of fault to x=0, dip-180 dip and negative sign on slip for surface dipping down to left
x5 = x5+xtip_axial_lower;

axial_lower_ends = [xtip_axial_lower 0; xtip_lower ztip_lower ];


%bottom axial shear zone
dip_axial_bot =  90 - (0 + dip_lower)/2;
slip = 2*srate*sin(abs(dip_lower)/2*pi/180);
xtip_axial_bot = xtip_lower+L_lower*cos(dip_lower*pi/180) + H/tan(dip_axial_lower*pi/180); 
L_axial_bot = H/sin(dip_axial_lower*pi/180);
[ux6,uz6,ux6e,uz6e,x6] = drive_visc_ramp_body_return_elastic(H,0,0,180-dip_axial_bot,L_axial_bot,-slip,mu,t,zobs); %shift tip of fault to x=0, dip-180 dip and negative sign on slip for surface dipping down to left
x6 = x6+xtip_axial_bot;

axial_bot_ends = [xtip_axial_bot 0; lower_ends(2,:)];




%interplate to xobs and add
[Xobs,Zobs] = meshgrid(xobs,zobs);
[X1,Z1] = meshgrid(x1,zobs);
[X2,Z2] = meshgrid(x2,zobs);
[X3,Z3] = meshgrid(x3,zobs);
[X4,Z4] = meshgrid(x4,zobs);
[X5,Z5] = meshgrid(x5,zobs);
[X6,Z6] = meshgrid(x6,zobs);

ux1 = griddata(X1,Z1,ux1,Xobs,Zobs);%,'nearest');
uz1 = griddata(X1,Z1,uz1,Xobs,Zobs);%,'nearest');
ux2 = griddata(X2,Z2,ux2,Xobs,Zobs);%,'nearest');
uz2 = griddata(X2,Z2,uz2,Xobs,Zobs);%,'nearest');
ux3 = griddata(X3,Z3,ux3,Xobs,Zobs);%,'nearest');
uz3 = griddata(X3,Z3,uz3,Xobs,Zobs);%,'nearest');
ux4 = griddata(X4,Z4,ux4,Xobs,Zobs);%,'nearest');
uz4 = griddata(X4,Z4,uz4,Xobs,Zobs);%,'nearest');
ux5 = griddata(X5,Z5,ux5,Xobs,Zobs);%,'nearest');
uz5 = griddata(X5,Z5,uz5,Xobs,Zobs);%,'nearest');
ux6 = griddata(X6,Z6,ux6,Xobs,Zobs);%,'nearest');
uz6 = griddata(X6,Z6,uz6,Xobs,Zobs);%,'nearest');

ux1e = griddata(X1,Z1,ux1e,Xobs,Zobs);%,'nearest');
uz1e = griddata(X1,Z1,uz1e,Xobs,Zobs);%,'nearest');
ux2e = griddata(X2,Z2,ux2e,Xobs,Zobs);%,'nearest');
uz2e = griddata(X2,Z2,uz2e,Xobs,Zobs);%,'nearest');
ux3e = griddata(X3,Z3,ux3e,Xobs,Zobs);%,'nearest');
uz3e = griddata(X3,Z3,uz3e,Xobs,Zobs);%,'nearest');
ux4e = griddata(X4,Z4,ux4e,Xobs,Zobs);%,'nearest');
uz4e = griddata(X4,Z4,uz4e,Xobs,Zobs);%,'nearest');
ux5e = griddata(X5,Z5,ux5e,Xobs,Zobs);%,'nearest');
uz5e = griddata(X5,Z5,uz5e,Xobs,Zobs);%,'nearest');
ux6e = griddata(X6,Z6,ux6e,Xobs,Zobs);%,'nearest');
uz6e = griddata(X6,Z6,uz6e,Xobs,Zobs);%,'nearest');


%% deep flat dislocation for elastic model
m=[lower_ends(2,1) -H 0 srate*t];
[ux7e, uz7e] =EdgeDisp(m,Xobs,-Zobs,0.25);



ux = ux1 + ux2 + ux3 + 0*(ux4 + ux5 + ux6);
uz = uz1 + uz2 + uz3 + 0*(uz4 + uz5 + uz6);

uxe = ux1e + ux2e + ux3e + 0*(ux4e + ux5e + ux6e) + ux7e;
uze = uz1e + uz2e + uz3e + 0*(uz4e + uz5e + uz6e) + uz7e;

% ux = ux1 + ux2 + ux3 + ux4 + ux5 + ux6;
% uz = uz1 + uz2 + uz3 + uz4 + uz5 + uz6;
% 
% uxe = ux1e + ux2e + ux3e + ux4e + ux5e + ux6e + ux7e;
% uze = uz1e + uz2e + uz3e + uz4e + uz5e + uz6e + uz7e;

figure
quiver(Xobs,-Zobs,ux,uz)
hold on
plot(upper_ends(:,1),-upper_ends(:,2),'r')
plot(middle_ends(:,1),-middle_ends(:,2),'r')
plot(lower_ends(:,1),-lower_ends(:,2),'r')
plot(axial_lower_ends(:,1),-axial_lower_ends(:,2),'k')
plot(axial_up_ends(:,1),-axial_up_ends(:,2),'k')
plot(axial_bot_ends(:,1),-axial_bot_ends(:,2),'k')
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
plot(axial_lower_ends(:,1),-axial_lower_ends(:,2),'k')
plot(axial_up_ends(:,1),-axial_up_ends(:,2),'k')
plot(axial_bot_ends(:,1),-axial_bot_ends(:,2),'k')
axis equal


figure
subplot(211)
plot(xobs,uxe(1,:))
grid on
subplot(212)
plot(xobs,uze(1,:))
grid on
