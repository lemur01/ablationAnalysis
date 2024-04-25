%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find all pixels on a line given by two points
% L. Muresan, lam94
% Bresenham's line algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function id = LineFind(x0, y0, x1, y1) 

steep = abs(y1 - y0) > abs(x1 - x0);
if steep 
     aux = x0; x0 = y0; y0 = aux;
     aux = x1; x1 = y1; y1 = aux;     
end
if x0 > x1 
    aux = x0; x0 = x1; x1 = aux;
    aux = y0; y0 = y1; y1 = aux;     
end
deltax = (x1 - x0);
deltay = abs(y1 - y0);
error = deltax / 2;

y = y0;
id = [];
if y0 < y1 ystep = 1; 
else ystep = -1;
end

for x = x0:x1
     if steep id = [id; y,x]; else id = [id; x,y]; end
     error = error - deltay;
     if error < 0
         y = y + ystep;
         error = error + deltax;
     end
     end
end