function pipesurface(portnumber)
%
%  Waits at socket portnumber ( between 5000 and 5010 ) until a
%  matrix arrives and then displays it in a surf.
% 
figure(1);
clg;
pause;
if ( nargin == 0 ) portnumber = 5001; end;
for i=0:10000,
  m = receive(portnumber);
  if ( sum(size(m)) ~= 0 ) 
     figure (1); plot(m(1,:)); pause(0.25);
  else
     break;
  end;
end;
