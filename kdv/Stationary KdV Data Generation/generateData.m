function generateData(func,xMin,xMax,npoints,fname)
%To generate waves, use cnoidalData or solitonData: this is a helper
%function
xrange = xMax-xMin;
dx = xrange / (npoints -1);
xs = xMin:dx:xMax;
vals = func(xs);
data = [xs;vals];
plot(xs,vals)
pause
fileID = fopen(fname,'w');
fprintf(fileID,'%6s %12s\n','x','f(x)');
fprintf(fileID,'%6.16f %12.16f\n',data);
fclose(fileID);
end