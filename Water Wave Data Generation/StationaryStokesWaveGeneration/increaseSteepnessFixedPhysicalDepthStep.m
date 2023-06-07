function [outname] = increaseSteepnessFixedPhysicalDepthStep(directory,fname,HTarget,Htol,cstep,specTol,cgTol,LyTol) 
fullFname = directory + fname;
[c1,N_org,h,~,y] = loadStokesSmall(fullFname);
C = c1 + cstep;
initRealH = realH(y,h);
initDeltaH = initRealH-HTarget;
newh = h;
[u,x_tilde,y1,HL,c,~,h,outSpec] = stokes_cg(y, N_org, C, newh,cgTol,LyTol);
if outSpec > specTol
    fprintf("Increasing number of points\n");
    N = 2 * N_org;
    N_sml = N_org;
    yn = y;
	yk = zeros(N, 1);
	ynk = fft(yn)/N_sml;
	yk(1:N_sml/2) = ynk(1:N_sml/2);
    %size(yk(3*N_sml/2+1:N));
    %size(ynk(N_sml/2+1:N_sml));
	yk(3*N_sml/2+1:N) = ynk(N_sml/2+1:N_sml);
	y = real(ifft(yk)*N);
    [u,x_tilde,y1,HL,~,~,h,~] = stokes_cg(y, N, C, newh,cgTol,LyTol);
    newRealH1 = realH(y1,h);
    deltaH1 = newRealH1-HTarget;
    if abs(deltaH1) < Htol
        %Save data
        hOpt = h;
        outH = realH(y,hOpt);
        outname = sprintf("H%3.10fHL%3.16f.txt",outH, HL/(2*pi));
        write_stokes(directory + outname, u', x_tilde, y, HL, C, N, h)
    else
        %Find a better h
        [smallh,smallDH,bigh,bigDH] = getToBisection(y,y1,N,C,h,HTarget,cgTol,LyTol);
        minDH = min(smallDH,bigDH);
        %Run bisection
        while abs(minDH) > Htol
            avgh = (bigh + smallh)/2;
            [~,~,ytest,~,~,~,~,~] = stokes_cg(y, N, C, avgh,cgTol,LyTol);
            testRealH = realH(ytest,avgh);
            testDh = testRealH - HTarget;
            if 0 > testDh 
                smallh = avgh;
                smallDH = testDh;
            else
                bigh = avgh;
                bigDH = testDh;
            end
            minDH = min(smallDH,bigDH)
            maxDH = max(smallDH,bigDH)
        end
        if minDH == smallDH
            %save small stokes
            hOpt = smallh;
            [u,x_tilde,yout,HL,c,~,h,~] = stokes_cg(y, 2 * N_org, C, hOpt,cgTol,LyTol);
            outH = realH(yout,h);
            outname = sprintf("H%3.10fHL%3.16f.txt",outH, HL/(2*pi));
            write_stokes(directory + outname, u', x_tilde, yout, HL, c, 2 * N_org, h)
        else
            %save big stokes
            hOpt = largeh;
            [u,x_tilde,yout,HL,c,~,h,~] = stokes_cg(y, 2 * N_org, C, hOpt,cgTol,LyTol);
            outH = realH(yout,h);
            outname = sprintf("H%3.10fHL%3.16f.txt",outH, HL/(2*pi));
            write_stokes(directory + outname, u', x_tilde, yout, HL, c, 2 * N_org, h)
        end
    end
else
    newRealH1 = realH(y1,h);
    deltaH1 = newRealH1-HTarget;
    if abs(deltaH1) < Htol
        %Save data
        hOpt = h;
        %[u,x_tilde,y,HL,c,~,h,~] = stokes_cg(directory, fname, N_org, C, hOpt,cgTol,LyTol);
        outH = realH(y,hOpt);
        outname = sprintf("H%3.10fHL%3.16f.txt",outH, HL/(2*pi));
        write_stokes(directory + outname, u', x_tilde, y1, HL, c, N_org, h)
    else
        %Find a better h
        [smallh,smallDH,bigh,bigDH] = getToBisection(y,y1,N_org,C,h,HTarget,cgTol,LyTol);
        minDH = min(smallDH,bigDH)
        maxDH = max(smallDH,bigDH)
        %Run bisection
        while abs(minDH) > Htol
            avgh = (bigh + smallh)/2;
            [~,~,ytest,~,~,~,~,~] = stokes_cg(y, N_org, C, avgh,cgTol,LyTol);
            testRealH = realH(ytest,avgh);
            testDh = testRealH - HTarget;
            if 0 > testDh 
                smallh = avgh;
                smallDH = testDh;
            else
                bigh = avgh;
                bigDH = testDh;
            end
            minDH = min(smallDH,bigDH)
            maxDH = max(smallDH,bigDH)
        end
        if minDH == smallDH
            %save small stokes
            hOpt = smallh;
            [u,x_tilde,yout,HL,c,~,h,~] = stokes_cg(y, N_org, C, hOpt,cgTol,LyTol);
            outH = realH(yout,h);
            outname = sprintf("H%3.10fHL%3.16f.txt",outH, HL/(2*pi));
            write_stokes(directory + outname, u', x_tilde, yout, HL, c, N_org, h)
        else
            %save big stokes
            hOpt = largeh;
            [u,x_tilde,yout,HL,c,~,h,~] = stokes_cg(y, N_org, C, hOpt,cgTol,LyTol);
            outH = realH(yout,h);
            outname = sprintf("H%3.10fHL%3.16f.txt",outH, HL/(2*pi));
            write_stokes(directory + outname, u', x_tilde, yout, HL, c, N_org, h)
        end
    end
end
end

function [smallh,smallDH,bigh,bigDH] = getToBisection(y, y1,n,c,newh,HTarget,cgTol,LyTol)
    newRealH1 = realH(y1,newh);
    deltaH1 = newRealH1-HTarget;
    newh2 = newh - deltaH1;
    [~,~,y2,~,~,~,~,~] = stokes_cg(y, n, c, newh2,cgTol,LyTol);
    newRealH2 = realH(y2,newh2);
    deltaH2 = newRealH2-HTarget;
    if (deltaH1 > 0) && (deltaH2 < 0)
        bigh = newh;
        bigDH = deltaH1;
        smallh = newh2;
        smallDH = deltaH2;
        else
            if (deltaH1 < 0) && (deltaH2 > 0)
               smallh = newh;
               smallDH = deltaH1;
               bigh = newh2;
               bigDH = deltaH2;
               %H1.0000000000HL0.0940093947879626 -> H1.0000000000HL0.0946708602691384
               %H1.5000000000HL0.1167600338982493 -> H1.5000000000HL0.1172090352464439
               %H2.0000000000HL0.1256740663114490 -> H2.0000000000HL0.1259727590558005
            else
                if deltaH1 < 0
                    i = 2;
                    while deltaH2 < 0
                        newh2 = newh - i * deltaH1;
                        [~,~,y2,~,~,~,~,~] = stokes_cg(y, n, c, newh2,cgTol,LyTol);
                        newRealH2 = realH(y2,newh2);
                        deltaH2 = newRealH2-HTarget;
                        i = i + 1
                    end
                    smallh = newh;
                    smallDH = deltaH1;
                    bigh = newh2;
                    bigDH = deltaH2;
                else
                    i = 2;
                    while deltaH2 > 0
                        newh2 = newh - i * deltaH1;
                        [~,~,y2,~,~,~,~,~] = stokes_cg(y, n, c, newh2,cgTol,LyTol);
                        newRealH2 = realH(y2,newh2);
                        deltaH2 = newRealH2-HTarget;
                        i = i + 1
                    end
                    smallh = newh2;
                    smallDH = deltaH2;
                    bigh = newh;
                    bigDH = deltaH1;
                end
            end
    end
end

function [H] = realH(y,h)
yk = fft(y);
N = length(yk);
yhat = yk(1)/N;
H = h - yhat;
end

function [u,x_tilde,y,HL,c,N,h,outSpec] = stokes_cg(yn, numit, speed, H,cgTol,LyTol)
global g c h N N_org
maxNumCompThreads(1);
	N = numit; c = speed;
	%N = 1048576; c = 1.0470936; 
	N_sml = N_org; %original # of points
	fprintf("Stokes speed c = %12.8e\n", c); 
	h = H; g = 1; 
    
	tol = cgTol; % conj-grad tolerance
	TOL = LyTol; % Newton-CG tolerance (if ||L y|| < TOL converged)
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
