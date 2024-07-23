function [H,V,converged,iters,err,derr] = wilson(S,maxiters,errtol,verb)

% Alternative version of Wilson's spectral factorisation algorithm (see
% also cpsd_specfac).
%
% Based on original code by M. Dhamala (mdhamala@bme.ufl.edu) &
% G. Rangarajan (rangaraj@math.iisc.ernet.in), UF, Aug 3-4, 2006.
%
% Reference: G. T. Wilson, The Factorization of Matricial Spectral
% Densities, SIAM J. Appl. Math, Vol. 23(4), pp 420-426, 1972.

if nargin < 2 || isempty(maxiters), maxiters = 100;  end
if nargin < 3 || isempty(errtol),   errtol   = 1e-9; end
if nargin < 4 || isempty(verb),     verb     = 0;    end

[n,~,h] = size(S);
h2 = 2*(h-1);

SX = cat(3,S,conj(S(:,:,h-1:-1:2))); % extend CPSD

% initialise P0
gamma =  ifft(SX,[],3);
gamma0 = gamma(:,:,1);
gamma0 = real((gamma0+gamma0')/2); % enforce symmetry, zero out imaginaries on diagonal
P0 = chol(gamma0);

P = repmat(P0,[1,1,h]);
P = cat(3,P,conj(P(:,:,h-1:-1:2)));

I = eye(n);

U = zeros(size(SX));
for i = 1:h2
    U(:,:,i) = chol(SX(:,:,i),'lower');
end

g = zeros(n,n,h2);
SF = zeros(n,n,h);

err = Inf;
for iters = 1:maxiters

    if verb, fprintf('iteration %3d ...',iters); end

    for i = 1:h2
        % Equivalent to: g(:,:,i) = P(:,:,i)\SX(:,:,i)/P(:,:,i)' + I;
        C = P(:,:,i)\U(:,:,i);
        g(:,:,i) = C*C' + I;
    end

    % []+ operation
    gam = real(ifft(g,[],3));
	gam(:,:,1) = gam(:,:,1)/2; % take half of the zero lag
	gp0 = gam(:,:,1);
	gam(:,:,h+1:end) = 0;      % zero out negative powers.
	gp = fft(gam,[],3);        % reconstitute

    T = -tril(gp0,-1);
    T = T-T';

    for i = 1:h2
        P(:,:,i) = P(:,:,i)*(gp(:,:,i) + T);
    end

    P0 = P0*(gp0 + T);

    % Calculate error
	for i = 1:h
		SF(:,:,i) = P(:,:,i)*P(:,:,i)';
	end
	olderr = err;
	err = max(abs(S(:)-SF(:)));
	derr = abs(err-olderr);

	% Check convergence
	if  err <= errtol
		converged = 1;
		if verb
			fprintf(' converged: err = %.2e, derr = %.2e\n',err,derr);
		end
		break
	end
	if  derr <= errtol
		converged = 2;
		if verb
			fprintf(2,' converged: err = %.2e, derr = %.2e - WARNING: tolerance not met\n',err,derr);
		else
			fprintf(2,'WARNING: converged in %d iterations, but tolerance not met: err = %.2e, derr = %.2e\n',iters,err,derr);
		end
		break
	end
	converged = 0;
	if verb
		fprintf(' err = %.2e, derr = %.2e\n',err,derr);
	end
end

if converged == 0
	fprintf(2,'WARNING - unconverged (exceeded maximum iterations)\n');
end

H = zeros(n,n,h);
for i = 1:h
    H(:,:,i) = P(:,:,i)/P0;
end

V = P0*P0';
