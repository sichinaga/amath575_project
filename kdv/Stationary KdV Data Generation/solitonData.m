function solitonData(a,t,xMin,xMax,npoints,fname)
%Generates data for soliton solution of KdV. a is a parameter, corresponding to amplitude and speed, t is time.
%xMin, xMax, and npoints determine the domain on which we are looking at
%our solution, and fname is the filename we'll be saved under.
    f = @(X) solitonFunc(a,X,t);
    generateData(f,xMin,xMax,npoints,fname)
end

function [out] = solitonFunc(a,x,t)
out = (4/3)*a^2*sech(a*(x-(2*a^2*t/3))).^2;
end