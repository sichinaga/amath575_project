function [fnames] = increaseDepthFixedSteepness(directory,fname,HLTarget,specTol,HLtol,cstep,outh,nSteps)
[~,~,h1,~,~] = loadStokesSmall(directory + fname);
hDiff = outh - h1;
hstep = hDiff/nSteps;
initGuessFileName = fname;
fnames = string(1:(nSteps-1));
for i = 1:nSteps
    initGuessFileName = increaseDepthFixedSteepnessStep(directory,initGuessFileName,HLTarget,HLtol,hstep,cstep,specTol);
    fnames(i) = initGuessFileName;
end
end

function [outname] = increaseDepthFixedSteepnessStep(directory,fname,HLTarget,HLtol,hstep,cstep,specTol) 
fullFname = directory + fname;
[c1,N_org,h1,~,~] = loadStokesSmall(fullFname);
newh = h1 + hstep;
[~,~,~,HL1,~,~,~,outSpec] = stokes_cg(directory, fname, N_org, c1, newh);
if outSpec > specTol
    [~,~,~,HL1,~,~,~,~] = stokes_cg(directory, fname, 2 * N_org, c1, newh);
    deltaHL1 = HL1-HLTarget;
    if abs(deltaHL1) < HLtol
        %Save data
        cOpt = c1;
        [u,x_tilde,y,HL,C,~,h,~] = stokes_cg(directory, fname, N_org, cOpt, newh);
        outH = realH(y,h);
        outname = sprintf("H%3.10fHL%3.16f.txt",outH, HL/(2*pi));
        write_stokes(directory + outname, u', x_tilde, y, HL, C, 2 * N_org, h)
    else
        %Find a better h
        [smallc,smallDHL,bigc,bigDHL] = getToBisection(directory, fname,2 * N_org,c1,newh,cstep,HLTarget);
        minDHL = min(smallDHL,bigDHL);
        %Run bisection
        while abs(minDHL) > HLtol
            avgc = (smallc + bigc)/2;
            [~,~,~,HLTest,~,~,~,~] = stokes_cg(directory, fname, 2 * N_org, avgc, newh);
            testDHL = HLTest - HLTarget;
            if 0 > testDHL 
                smallc = avgc;
                smallDHL = testDHL;
            else
                bigc = avgc;
                bigDHL = testDHL;
            end
            minDHL = min(smallDHL,bigDHL);
        end
        if minDHL == smallDHL
            %save small stokes
            cOpt = smallc;
            [u,x_tilde,y,HL,C,~,h,~] = stokes_cg(directory, fname, 2 * N_org,cOpt,newh);
            outH = realH(y,h);
            outname = sprintf("H%3.10fHL%3.16f.txt",outH, HL/(2*pi));
            write_stokes(directory + outname, u', x_tilde, y, HL, C, 2 * N_org, h)
        else
            %save big stokes
            cOpt = largec;
            [u,x_tilde,y,HL,C,~,h,~] = stokes_cg(directory, fname, 2 * N_org, cOpt, newh);
            outH = realH(y,h);
            outname = sprintf("H%3.10fHL%3.16f.txt",outH, HL/(2*pi));
            write_stokes(directory + outname, u', x_tilde, y, HL, C, 2 * N_org, h)
        end
    end
else
    deltaHL1 = HL1-HLTarget;
    if abs(deltaHL1) < HLtol
        %Save data
        cOpt = c1;
        [u,x_tilde,y,HL,C,~,h,~] = stokes_cg(directory, fname, N_org, cOpt, hNew);
        outH = realH(y,h);
        outname = sprintf("H%3.10fHL%3.16f.txt",outH, HL/(2*pi));
        write_stokes(directory + outname, u', x_tilde, y, HL, C, N_org, h)

    else
        %Find a better h
        [smallc,smallDHL,bigc,bigDHL] = getToBisection(directory, fname, N_org,c1,newh,cstep,HLTarget);
        minDHL = min(smallDHL,bigDHL);
        %Run bisection
        while abs(minDHL) > HLtol
            avgc = (smallc + bigc)/2;
            [~,~,~,HLTest,~,~,~,~] = stokes_cg(directory, fname, N_org, avgc, newh);
            testDHL = HLTest - HLTarget;
            if 0 > testDHL 
                smallc = avgc;
                smallDHL = testDHL;
            else
                bigc = avgc;
                bigDHL = testDHL;
            end
            minDHL = min(smallDHL,bigDHL);
        end
        if minDHL == smallDHL
            %save small stokes
            cOpt = smallc;
            [u,x_tilde,y,HL,C,~,h,~] = stokes_cg(directory, fname, N_org,cOpt,newh);
            outH = realH(y,h);
            outname = sprintf("H%3.10fHL%3.16f.txt",outH, HL/(2*pi));
            write_stokes(directory + outname, u', x_tilde, y, HL, C, N_org, h)
        else
            %save big stokes
            cOpt = largec;
            [u,x_tilde,y,HL,C,~,h,~] = stokes_cg(directory, fname, N_org, cOpt, newh);
            outH = realH(y,h);
            outname = sprintf("H%3.10fHL%3.16f.txt",outH, HL/(2*pi));
            write_stokes(directory + outname, u', x_tilde, y, HL, C, N_org, h)
        end
    end
end
end

function [smallc,smallDHL,bigc,bigDHL] = getToBisection(directory, fname, n,c1, newh,cstep,HLTarget)
    [~,~,~,HL1,~,~,~,~] = stokes_cg(directory, fname, n, c1, newh);
    deltaHL1 = HL1-HLTarget;
    if deltaH1 * cstep < 0
        cstep = -cstep;
    end
    c2 = c1 - cstep;
    [~,~,~,HL2,~,~,~,~] = stokes_cg(directory, fname, n, c2, newh);
    deltaHL2 = HL2-HLTarget;
    i = 2;
    while deltaHL2 * deltaHL1 > 0
        c2 = c1 - i * cstep;
        [~,~,~,HL2,~,~,~,~] = stokes_cg(directory, fname, n, c2, newh);
        deltaHL2 = HL2-HLTarget;
        i = i + 1;
    end

    if deltaHL1 > 0
        smallc = c2;
        smallDHL = deltaHL2;
        bigc = c1;
        bigDHL = deltaHL1;
    else
        smallc = c1;
        smallDHL = deltaHL1;
        bigc = c2;
        bigDHL = deltaHL2;
    end
end

function [H] = realH(y,h)
yk = fft(y);
N = length(yk);
yhat = yk(1)/N;
H = h - yhat;
end

function [u,x_tilde,y,HL,c,N,h,outSpec] = stokes_cg(directory, fname, numit, speed, H)
global g c h N N_org
maxNumCompThreads(1);
    fullFname = directory + fname;
	[yn, z_tilde_id, u_old] = load_stokes(fullFname);
	N = numit; c = speed;
	%N = 1048576; c = 1.0470936; 
	N_sml = N_org; %original # of points
	fprintf("Stokes speed c = %12.8e\n", c); 
	h = H; g = 1; 
	if (N ~= N_sml)
		fprintf("Increasing number of points\n");
		yk = zeros(N, 1);
		ynk = fft(yn)/N_sml;
		
		yk(1:N_sml/2) = ynk(1:N_sml/2);
        size(yk(3*N_sml/2+1:N))
        size(ynk(N_sml/2+1:N_sml))
		yk(3*N_sml/2+1:N) = ynk(N_sml/2+1:N_sml);
		yn = real(ifft(yk)*N);
	end
    
	tol = 1e-14; % conj-grad tolerance
	TOL = 1e-10; % Newton-CG tolerance (if ||L y|| < TOL converged)
	maxiter = 10000; % max # of CG iterations 
	l = 0;
	k = fftshift(-N/2:(N/2-1)); k = transpose(k);
	k_op = k.*coth(k*h);
        k_op(1) = 0; %Maybe k_op(1) = 1./h;	
	while 1
		b = -stokes_eq(yn, k_op);
		dy = conjgrad(b, 0*yn, yn, tol, k_op, maxiter);
		yn = yn+dy;
		yn = real(yn);
		l = l+1;
		fprintf("Ly = %23.16e at iter = %d\n", norm(stokes_eq(yn, k_op)), l);
		if (norm(stokes_eq(yn, k_op)))<=TOL
			break;
		end
	end
	% Checking Fourier spectrum
	spec = fft(yn);
    outSpec = spec(N/2+1);
	fprintf("Largest negative Fourier harmonic #%d\tyk = %12.8e\n", k(N/2+1), spec(N/2+1)); 
	% Computing new x_tilde, y
	y = yn;
	T = -1i*coth(k*h);
        T(1) = 0;	
	x_tilde = ifft(T.*fft(y)); 
	z_tilde = x_tilde+1i*y;
	q = pi*(2*(0:N-1)/N - 1);
	u = 2.*atan2(sin(0.5*q), cos(0.5*q));
	HL = max(y)-min(y);
end

function x = conjgrad(b, x, yn, tol, k_op, maxiter)
	r = b-matrix_mult(x, yn, k_op);
	rprodold = r'*r;
	if sqrt(rprodold) <= tol
		x = x;
        l=0;
	else
		p = r;
		l = 0;
		for m = 1:maxiter
			Ap = matrix_mult(p, yn, k_op);
			alpha = rprodold./(p'*Ap);
			x = x+alpha*p;
			r = r-alpha*Ap;
			rprodnew = r'*r;
			l = l+1;
			if sqrt(rprodnew) <= tol
				break;
			end
			beta = rprodnew./rprodold;
			p = r+beta*p;
			rprodold = rprodnew;
		end
	end
	fprintf("Residual = %23.16e at iter =  %d\n", sqrt(rprodold), l);
end

function Av = matrix_mult(v, yn, k_op)
global g c h N N_org
	vk = fft(v);
	yk = fft(yn);

	Av = (c^2/g)*ifft(k_op.*vk)-v-(ifft(k_op.*fft(yn.*v))+v.*ifft(k_op.*yk)+yn.*ifft(k_op.*vk));
end

function Ly = stokes_eq(v, k_op)
global g c h N N_org
	vk = fft(v);
	
	Ly = (c^2/g)*ifft(k_op.*vk)-v-(ifft(k_op.*fft(0.5*v.*v))+v.*ifft(k_op.*vk));
end

function [y, z_tilde, u_old, N_org] = load_stokes(stokes_file)
global g c h N N_org
  fhin = fopen(stokes_file, "r");
  if (fhin == -1)
    fprintf("Cannot open Stokes wave file\n");
    exit(0); 
  else 
    line = fgets(fhin);
    line = fgets(fhin);
    
    % Matlab
    v = sscanf(line, "%%# N = %d\tL = %f\n");    
    N_org = v(1);

    fclose(fhin); 
  end
  raw = load(stokes_file);
  z_tilde = raw(:,2) + 1.i*raw(:,3);
  y = imag(z_tilde);

  u_old = raw(:,1);
end

function [c,N,h,x,y] = loadStokesSmall(fname)
fh = fopen(fname,"r");
if (fh == -1)
    fprintf("Cannot read the stokes wave file\n");
    exit(0);
else
    line = fgets(fh);
    line = fgets(fh);
    vec = sscanf(line,"%%# N = %d\tL = %d\n");
    N = vec(1);
    line = fgets(fh);
    data = sscanf(line,"%%# H/L = %f\tc = %f\th = %f\ty0 = %f\tOmega = %f\n");
    c =  data(2);
    h = data(3);
    fclose(fh);
end
rawData = load(fname);
u = rawData(:,1);
xtilde = rawData(:,2);
x = u + xtilde;
y = rawData(:,3);
end


function write_stokes(fname, u, x_tilde, y, HL, c, N, h)
  fh = fopen(fname,"w");
  fprintf(fh, "%%# 1. u 2. x_tilde 3. y \n");
  fprintf(fh, "%%# N = %d\tL = %23.16e\n", N, 0.);
  fprintf(fh, "%%# H/L = %23.16e\tc = %23.16e\th = %23.16e\ty0 = %e\tOmega = %e\n\n", HL, c, h, 0., 0.);
  for j = 1:length(y)
    fprintf(fh, "%23.16e\t%23.16e\t%23.16e\n", u(j), x_tilde(j), y(j));
  end
  fclose(fh);
end
