	function [u1, u2] =EdgeDisp(m,x1,x2,nu)
%EDGEDISP    [u1, u2] =EdgeDisp(m,x1,x2,nu);
 
%Computes displacements at 'x' caused by the edge dislocation
%specified in 'm':
%     m(1) = Horizontal position of updip end
%     m(2) = Depth of updip end (MUST be negative)
%     m(3) = Dip (degrees) measured positve down from x1 axis;
%     m(4) = Slip 
%
%  Observation coordinates, 'x1', and 'x2' can be a matrix;
%  'nu' is Poisson's ratio, which if omitted defaults to 0.25.
%  Output is two matricies 'u1' and 'u2' same size as x1 and x2
%  in same Units as the dislocation
%
%
%February 10, 2000, P. Cervelli
%Revised April 17, P. Segall

%Check for omitted Poisson's ratio   
if nargin < 3      nu=0.25;   end

% components of slip vector
	b1=m(4)*cos(m(3)*pi/180);  
	b2=m(4)*sin(m(3)*pi/180); 

% Some parameters used in calcaulating displacements   
	C=1/(pi*(1-nu));   
	r1=sqrt((x2-m(2)).^2+(x1-m(1)).^2);   
	r2=sqrt((x2+m(2)).^2+(x1-m(1)).^2);   
	dx1=x1-m(1);   
	dx2=x2-m(2);   
	dx2a=x2+m(2); 
 
% define branch cuts
	% rotated coordinates at source
	x1p = cos(m(3)*pi/180)*dx1 - sin(m(3)*pi/180)*dx2;
	x2p = cos(m(3)*pi/180)*dx2 + sin(m(3)*pi/180)*dx1;

	% rotated coordinates at image
	x1pi = cos(m(3)*pi/180)*dx1 - sin(m(3)*pi/180)*dx2a;
	x2pi = cos(m(3)*pi/180)*dx2a + sin(m(3)*pi/180)*dx1;

	theta1=atan2(x2p,-x1p);   
	theta2=-atan2(x2pi,x1pi); 

  
u1 =  b1*C*( (1-nu)/2*(theta2-theta1) + dx1.*dx2./(4*r1.^2) - ...  
	   dx1.*( x2+(3-4*nu)*m(2) )./(4*r2.^2) + ...  
	   (x2.*m(2).*dx1.*dx2a)./r2.^4 ) + ...  
      b2*C*( (1-2*nu)/4*(log(r2) - log(r1)) - dx2.^2./(4*r1.^2) + ...  
	  (x2.^2 + 4*(nu-1)*x2*m(2) + (4*nu-3)*m(2)^2 )./(4*r2.^2) + ... 	
	  (x2*m(2).*dx2a.^2)./r2.^4 );

u2 = b2*C*( (1-nu)/2*(theta1-theta2) + dx1.*dx2./(4*r1.^2) - ...    	
	   dx1.*( x2+(3-4*nu)*m(2))./(4*r2.^2) -...
	   dx1.*x2*m(2).*dx2a./r2.^4 ) + ...  
     b1*C*( (1-2*nu)/4*(log(r2) - log(r1)) +  dx2.^2./(4*r1.^2) - ...  
	    (dx2a.^2 - 2*m(2)^2 - 2*(1-2*nu)*m(2)*dx2a)./(4*r2.^2) + ...  
	    (x2*m(2).*dx2a.^2)./r2.^4 );
