function [fnames] = increaseSteepnessFixedPhysicalDepth(directory,fname,HTarget,specTol,Htol,outC,nSteps,cgTol,LyTol)
[c1,~,~,~,~] = loadStokesSmall(directory + fname);
cDiff = outC - c1;
deltaC = cDiff/nSteps;
initGuessFileName = fname;
fnames = string(1:(nSteps-1));
for i = 1:nSteps
    initGuessFileName = increaseSteepnessFixedPhysicalDepthStep(directory,initGuessFileName,HTarget,Htol,deltaC,specTol,cgTol,LyTol); 
    fnames(i) = initGuessFileName;
end
end
%Try depth 2
%Try fixed amplitude, varying Kh max, 
%Get oscillation
%Get corner angles