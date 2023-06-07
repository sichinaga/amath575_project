function cnoidalData(k,gamma,b,t,xMin,xMax,npoints,fname)
%Generates data for cnoidal solution of KdV. k, gamma, and b is are parameters,
% t is time. Requires symbolic computation package
%xMin, xMax, and npoints determine the domain on which we are looking at
%our solution, and fname is the filename we'll be saved under.
    f = @(X) cnoidal(k,gamma,b,t,X);
    generateData(f,xMin,xMax,npoints,fname)
end

function [out] = cnoidal(k,gamma,b,t,x)
a = 2*(k^2)*(gamma^2);
V =  6*b + 4*(2*(k^2)-1)*gamma;
out = b + a * jacobiCN(gamma*(x-V*t),k^2).^2;
end
