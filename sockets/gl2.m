function gl2(portnumber)
%
%  Waits at socket portnumber (between 5000 and 5010) until a matrix arrives
% 
%  This file is for the Ginzburg-Landau 2-D (gl2) problem
%

figure(1);
clg;
if ( nargin == 0 ) portnumber = 5001; end;
for i=0:10000
    m = receive(portnumber);
    if ( sum(size(m)) ~= 0 ) 
       n = size(m,1)
       nx = sqrt(n/4)
       ny = nx 
       for j = 1: nx
           for i = 1:ny
               v1(i,j) = m((i-1)*nx+j,1) ;
               v2(i,j) = m(nx*ny+(i-1)*nx+j,1) ;
           end
       end

%      Compute the magnitude

       for j = 1: nx
           for i = 1:ny
               v(i,j) = norm([v1(i,j),v2(i,j)]) ;
           end
       end

%      Surface plots

       zz = [ v v ; v v ];
       figure (1); surf(zz); pause (0.25);
       figure (2); contour(zz); pause (0.25);

   else
      break;
   end;
end;
