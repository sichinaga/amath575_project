function [outZ,outPsi] = WaterWaveForwards(initYs,nT,deltaT,c1,c2)
	%Uses RK4 to iterate the water wave problem forward in time given an initial surface elevation of two stokes waves set next to each other. Does not work
    h = 1;
    N = length(initYs);
    k = fftshift(-N/2:1:(N/2)-1)';
    That = -1i * coth(h * k);
    That(1) = 0;
    Rhat = 1i * tanh(h * k);
    du = 1i * k;
    g = 1;
    if (c1 == c2)
        c = c1;
        [startZ,startPsi]=startVals(initYs,That,c);
    else
        x = ifft(That .* fft(initYs));
        startZ = (x - 1i * initYs)';
        nHalf = N/2;
        Ksplit = fftshift(nHalf/2:1:(nHalf/2)-1)';
        ThatSplit = -1i * coth(h * Ksplit);
        ThatSplit(1) = 0;
        y1 = initYs(1:nHalf);
        y2 = initYs(nHalf+1:N);
        psi1 = c1 * ifft(ThatSplit .* fft(y1));
        psi2 = c2 * ifft(ThatSplit .* fft(y2));
        startPsi = [psi1 ; psi2];
    end
    outZ = zeros(length(startZ),nT+1);
    outPsi = zeros(length(startPsi),nT+1);
    outZ(:,1) = startZ;
    outPsi(:,1) = startPsi;
    for i = 1:1:(nT)
        [newZ,newPsi] = iterateStep(outZ(:,i),outPsi(:,i),That,Rhat,du,g,deltaT);
        outZ(:,i+1)=newZ;
        outPsi(:,i+1)=newPsi;
    end
end

function [startZ,startPsi] = startVals(inputY,That,c)
    x = ifft(That .* fft(inputY));
    startZ = x + 1i * inputY;
    startPsi = c * x;
    %is it the case that duTduf = Tduduf = duduTf
end

function [outZ,outPsi] = iterateStep(inZ,inPsi,That,Rhat,du,g,deltaT)
%Uses RK4
    H = deltaT/2;
    [z1,psi1] = derivs(inZ,inPsi,That,Rhat,du,g);
    [z2,psi2] = derivs(inZ + H*z1,inPsi + H*psi1,That,Rhat,du,g);
    [z3,psi3] = derivs(inZ + H*z2,inPsi + H*psi2,That,Rhat,du,g);
    [z4,psi4] = derivs(inZ + deltaT*z3,inPsi + deltaT*psi3,That,Rhat,du,g);
    outZ = inZ + (deltaT/6) * (z1 + z2 + z3 + z4);
    outPsi = inPsi + (deltaT/6) * (psi1 + psi2 + psi3 + psi4);  
end

function [zt,psit] = derivs(inZ,inPsi,That,Rhat,du,g)
    zu = ifft(du .* fft(inZ));
    normzu2 = real(zu).^2 + imag(zu).^2;
    psik = fft(inPsi);
    Psiu = ifft(du .* psik);
    rhatDuPsi = ifft(du .* Rhat .* psik);
    rhatDuPsidiv = rhatDuPsi ./ normzu2;
    ThatRduPsi = ifft(That .* fft(rhatDuPsidiv));
    zt = -zu .*(1i .* rhatDuPsidiv + ThatRduPsi);
    psiuThatpsiu = Psiu .* rhatDuPsi;
    psit1 = -(ifft(That .* psiuThatpsiu))./normzu2;
    psit2 = -Psiu .* ThatRduPsi;
    psit = -g*imag(inZ) + psit1 + psit2;
end