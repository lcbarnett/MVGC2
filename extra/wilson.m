function [H,V,converged,iters,aerr,derr] = wilson(S,maxi,tol,verb)

% Alternative version of Wilson's spectral factorisation algorithm (see
% also cpsd_specfac).
%
% Based on original code by M. Dhamala (mdhamala@bme.ufl.edu) &
% G. Rangarajan (rangaraj@math.iisc.ernet.in), UF, Aug 3-4, 2006.
%
% Reference: G. T. Wilson, The Factorization of Matricial Spectral
% Densities, SIAM J. Appl. Math, Vol. 23(4), pp 420-426, 1972.

if nargin < 2 || isempty(maxi), maxi = 100;   end
if nargin < 3 || isempty(tol),  tol  = 1e-10; end
if nargin < 4 || isempty(verb), verb = false; end

[n,~,h] = size(S);
h2 = 2*(h-1);

MAS = max(abs(S(:)));

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

converged = false;
aerr = Inf;
for iters = 1:maxi

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
    oerr = aerr;
	for i = 1:h
		SF(:,:,i) = P(:,:,i)*P(:,:,i)';
	end
	aerr = max(abs(SF(:)-S (:)))/MAS; % absolute error
	derr = abs(oerr-aerr);            % error difference

	% Check convergence
	if verb
		fprintf(' aerr = %.2e, derr = %.2e',aerr,derr);
	end
	if  derr <= tol
		converged = true;
		if verb
			fprintf(' - converged\n');
		end
		break
	end
	if verb
		if iters == maxi
			fprintf(' - unconverged (timed out)\n');
		else
			fprintf('\n');
		end
	end
end

H = zeros(n,n,h);
for i = 1:h
    H(:,:,i) = P(:,:,i)/P0;
end

V = P0*P0';
