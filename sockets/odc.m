function odc(portnumber)
%
%  Waits at socket portnumber (between 5000 and 5010) until a matrix arrives
% 
%  This file is for the optimal design with composites (odc) problem
%

figure(1);
clg;
if ( nargin == 0 ) portnumber = 5001; end;
for i=0:10000
    m = receive(portnumber);
    if ( sum(size(m)) ~= 0 ) 
       n = size(m,1)
       nx = sqrt(n)
       ny = n/nx 
       for j = 1: nx
           for i = 1:ny
               v(i+1,j+1) = m((i-1)*nx+j,1) ;
           end
       end
       for j = 1:nx+2
           v(1,j) = 0.0 ;
           v(ny+2,j) = 0.0 ;
       end
       for i = 1:ny+2
           v(i,1) = 0.0 ;
           v(i,nx+2) = 0.0 ;
       end

%      Surface plots

       hx = 1/(nx+1); hy = 1/(ny+1);
       [px,py] = gradient(v,hx,hy);
       r = sqrt(px.^2 + py.^2) ;
       figure(1); surf(r);  pause (0.25);
       figure (2); contour(r); pause(0.25);
   else
      break;
   end;
end;
