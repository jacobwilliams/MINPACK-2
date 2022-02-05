function fdc(portnumber)
%
%  Waits at socket portnumber (between 5000 and 5010) until a matrix arrives
% 
%  This file is for the flow in a channel (fdc) problem
%

figure(1);
clg;
figure(2);
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
               psi(i+1,j+1) = m((i-1)*nx+j,1) ;
           end
       end
       for j = 1:nx+2
           psi(1,j) = 0.0 ;
           psi(ny+2,j) = 0.0 ;
       end
       for i = 1:ny+2
           psi(i,1) = 0.0 ;
           psi(i,nx+2) = 0.0 ;
       end

%      Streamlines plots

       V = [ -.11 -.1 -.08 -.06 -.04 -.02 -.01 ...
            -1.0e-5 1.0e-6 1.0e-5 1.0e-4 1.0e-3 2.0e-3 ] ;
       disp(sprintf('Beginning contour plot'))
       figure(1); contour(psi,V);  pause (0.0);
       disp(sprintf('Finished contour plot'))

%      Equivorticity plots

       hx = 1/(nx+1); hy = 1/(ny+1);
       for j = 1:nx
           for i = 1:ny
               v(i,j) = -(psi(i+2,j+1)-2*psi(i+1,j+1)+psi(i,j+1))/(hx*hx) ...
                        -(psi(i+1,j+2)-2*psi(i+1,j+1)+psi(i+1,j))/(hy*hy) ;
           end
       end
       % V = [ -5.0 -3.0 -1.0 1.0 3.0 5.0 ] ;
       V = [ -5.0 -4.0 -3.0 -2.0 -1.0 0.0 1.0 2.0 3.0 4.0 5.0 ] ;
       disp(sprintf('Beginning contour plot'))
       figure(2); contour(v,V);  pause (0.0);
       disp(sprintf('Finished contour plot'))
   else
      break;
   end;
end;
