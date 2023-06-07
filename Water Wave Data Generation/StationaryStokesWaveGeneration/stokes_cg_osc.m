function [outname] = stokes_cg_osc(directory, fname1,fname2,cgtol,LyTol,maxiter)
global g c h N N_org
%fname1 is for closer to max, fname2 is right next to it
maxNumCompThreads(1);
    fullFname1 = directory + fname1;
    fullFname2 = directory + fname2;

	[y1, z_tilde_id1, u_old1] = load_stokes(fullFname1);
	[y2, z_tilde_id, u_old] = load_stokes(fullFname2);
    [C,~,H,~,~] = loadStokesSmall(fullFname1);
    sml_yn = y1-y2;
	yn = y1+sml_yn;
	N = 2 * N_org; c = C;
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
    
	tol = cgtol; % conj-grad tolerance
	TOL = LyTol; % Newton-CG tolerance (if ||L y|| < TOL converged)
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
	fprintf("Largest negative Fourier harmonic #%d\tyk = %12.8e\n", k(N/2+1), spec(N/2+1)); 
	if (spec(N/2+1)>1e-9)
		fprintf("Spectrum unresolved\n");
		exit(0); 
	end
	% Computing new x_tilde, y
	y = yn;
	T = -1i*coth(k*h);
        T(1) = 0;	
	x_tilde = ifft(T.*fft(y)); 
	z_tilde = x_tilde+1i*y;
	q = pi*(2*(0:N-1)/N - 1);
	u = 2.*atan2(sin(0.5*q), cos(0.5*q));
	HL = max(y)-min(y);
    outname = sprintf("h%3.16fHL%3.16f.txt", h, HL/(2*pi));
    filenameanddir = sprintf(directory+"h%3.16fHL%3.16f.txt", h, HL/(2*pi));
	write_stokes(filenameanddir, u', x_tilde, y, HL/(2*pi), c, N, h);
	%figure(1)
	%subplot(2,1,1)
	%plot(u'+x_tilde, y, u_old+real(z_tilde_id), imag(z_tilde_id))
	%subplot(2,1,2)
	%semilogy((-N/2:N/2-1), fftshift(abs(fft(y)))/N)
end

function x = conjgrad(b, x, yn, tol, k_op, maxiter)
	r = b-matrix_mult(x, yn, k_op);
	rprodold = r'*r;
	%rprodold
	%pause
	if sqrt(rprodold) <= tol
		x = x;
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
			%sqrt(rprodnew)
			%pause
			if sqrt(rprodnew) <= tol
				break;
			end
			beta = rprodnew./rprodold;
			p = r+beta*p;
			rprodold = rprodnew;
		end
	end
	fprintf("Residual = %23.16e at iter =  %d\n", sqrt(rprodnew), l);
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

function [y, z_tilde, u_old] = load_stokes(stokes_file)
global g c h N N_org
  fhin = fopen(stokes_file, "r");
  if (fhin == -1)
    fprintf("Cannot open Stokes wave file\n");
    exit(0); 
  else 
    line = fgets(fhin);
    line = fgets(fhin);
    %% Octave
    %[v1, v2] = sscanf(line, "%# N = %s\tL = %s\n","C");
    %N_org = str2num(v1); 
    %l  = str2double(v2);
    
    % Matlab
    v = sscanf(line, "%%# N = %d\tL = %f\n");    
    N_org = v(1);
    %l = v(2);

    %line = fgets(fhin);
    %% Octave
    %[v1, v2, v3, v4] = sscanf(line, "# H/L = %s\tc = %s\ty0 = %s\tOmega = %s\n","C");
    %c = str2num(v2);
    
    %fprintf("Stokes speed c = %12.8e\n", c); 
    fclose(fhin); 
  end
  raw = load(stokes_file);
  z_tilde = raw(:,2) + 1.i*raw(:,3);
  y = imag(z_tilde);

  %q = pi*(2*(0:N-1)/N - 1);
  %u = 2.*atan2(sin(0.5*q), cos(0.5*q));
  u_old = raw(:,1);
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

function check_doubling_points()
	%q = pi*(2*(0:N_sml-1)/N_,sml - 1);
	%u1 = 2.*atan2(sin(0.5*q), cos(0.5*q));
	%k = fftshift(-N/2:(N/2-1)); k = transpose(k);
	%T = -1i*coth(k*h);
        %T(1) = 0;	
	%x_tilde = ifft(T.*fft(yn)); 
	%q = pi*(2*(0:N-1)/N - 1);
	%u2 = 2.*atan2(sin(0.5*q), cos(0.5*q));
	%plot(u2'+x_tilde, yn, '+', u1'+x_tildeold, yold, 'o')
	%pause
end
