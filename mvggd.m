function y = mvggd(X, Mu, Sigma, Beta)
%MVGGD Multivariate Generalized Gaussian Distribution (GGD).
%   Y = MVGGD(X) returns the probability density of the multivariate
%   generalized Gaussian distribution with zero mean, identity covariance 
%   matrix, and unitary shaping, evaluated at each row of X.  Rows of the 
%   N-by-D matrix X correspond to observations or points, and columns 
%   correspond to variables or coordinates.  Y is an N-by-1 vector.
%
%   Y = MVGGD(X,MU) returns the density of the multivariate GGD with mean 
%   MU, identity covariance matrix, and BETA equal to 1, evaluated at each
%   row of X.  MU is a 1-by-D vector, or an N-by-D matrix, in which case 
%   the density is evaluated for each row of X with the corresponding row 
%   of MU.  MU can also be a scalar value, which MVGGD replicates to match 
%   the size of X.
%
%   Y = MVGGD(X,MU,SIGMA) returns the density of the multivariate GGD with 
%   mean MU, covariance SIGMA, and unitary shaping, evaluated at each row
%   of X.  SIGMA is a D-by-D matrix, or an D-by-D-by-N array, in which case
%   the density is evaluated for each row of X with the corresponding page
%   of SIGMA, i.e., MVGGD computes Y(I) using X(I,:) and SIGMA(:,:,I).
%   If the covariance matrix is diagonal, containing variances along the 
%   diagonal and zero covariances off the diagonal, SIGMA may also be
%   specified as a 1-by-D matrix or a 1-by-D-by-N array, containing 
%   just the diagonal. Pass in the empty matrix for MU to use its default 
%   value when you want to only specify SIGMA.
%
%   Y = MVGGD(X,MU,SIGMA,BETA) returns the density of the multivariate GGD 
%   with mean MU, covariance SIGMA, and BETA, evaluated at each row of X.  
%   BETA is a scalar, or an N-by-1 vector, in which case 
%   the density is evaluated for each row of X with the corresponding row 
%   of BETA.
%
%   If X is a 1-by-D vector, MVGGD replicates it to match the leading
%   dimension of MU or the trailing dimension of SIGMA. If BETA = 1, MVGGD
%   becomes MVNPDF
%
%   Example:
%
%      mu = [1 -1]; Sigma = [.9 .4; .4 .3]; Beta = 2.2; 
%      [X1,X2] = meshgrid(linspace(-2,5,25)', linspace(-3,1,25)');
%      X = [X1(:) X2(:)];
%      y = mvggd(X, mu, Sigma, Beta);
%      surf(X1,X2,reshape(y,25,25));


if nargin<1
    error(message('stats:mvggd:TooFewInputs'));
elseif ~ismatrix(X)
    error(message('stats:mvggd:InvalidData'));
end

% Get size of data.  Column vectors provisionally interpreted as multiple scalar data.
[n,d] = size(X);
if d<1
    error(message('stats:mvggd:TooFewDimensions'));
end

% Assume zero mean, data are already centered
if nargin < 2 || isempty(Mu)
    X0 = X;
    Sigma = eye(d); 

% Get scalar mean, and use it to center data
elseif isscalar(Mu)
    X0 = X - Mu;

% Get vector mean, and use it to center data
elseif ismatrix(Mu)
    [n2,d2] = size(Mu);
    if d2 ~= d % has to have same number of coords as X
        error(message('stats:mvggd:ColSizeMismatch'));
    elseif (n2 == 1) || (n2 == n) % mean is a single row or a full vector.
        X0 = X - Mu;
    elseif n == 1 % data is a single row, rep it out to match mean
        n = n2;
        X0 = X - Mu;  
    else % sizes don't match
        error(message('stats:mvggd:RowSizeMismatch'));
    end
    
else
    error(message('stats:mvggd:BadMu'));
end

% Assume identity covariance, data are already standardized
if nargin < 3 || isempty(Sigma)
    % Special case: if Sigma isn't supplied, then interpret X
    % and Mu as row vectors if they were both column vectors
    if (d == 1) && (numel(X) > 1)
        X0 = X0';
        d = size(X0,2);
    end
    xRinv = X0;
    Sigma = eye(d); 
    
% Single covariance matrix
elseif ismatrix(Sigma)
    sz = size(Sigma);
    if sz(1)==1 && sz(2)>1
        % Just the diagonal of Sigma has been passed in.
        sz(1) = sz(2);
        sigmaIsDiag = true;
    else
        sigmaIsDiag = false;
    end
    
    % Special case: if Sigma is supplied, then use it to try to interpret
    % X and Mu as row vectors if they were both column vectors.
    if (d == 1) && (numel(X) > 1) && (sz(1) == n)
        X0 = X0';
        d = size(X0,2);
    end
    
    %Check that sigma is the right size
    if sz(1) ~= sz(2)
        error(message('stats:mvggd:BadCovariance'));
    elseif ~isequal(sz, [d d])
        error(message('stats:mvggd:CovSizeMismatch'));
    else
        if sigmaIsDiag
            if any(Sigma<=0)
                error(message('stats:mvggd:BadDiagSigma'));
            end
            R = sqrt(Sigma);
            xRinv = X0./R;
        else
            % Make sure Sigma is a valid covariance matrix
            [R,err] = cholcov(Sigma,0);
            if err ~= 0
                error(message('stats:mvggd:BadMatrixSigma'));
            end
            % Create array of standardized data, and compute log(sqrt(det(Sigma)))
            xRinv = X0 / R;
        end
    end
    
% Multiple covariance matrices
elseif ndims(Sigma) == 3
    
    sz = size(Sigma);
    if sz(1)==1 && sz(2)>1
        % Just the diagonal of Sigma has been passed in.
        sz(1) = sz(2);
        Sigma = reshape(Sigma,sz(2),sz(3))';
        sigmaIsDiag = true;
    else
        sigmaIsDiag = false;
    end

    % Special case: if Sigma is supplied, then use it to try to interpret
    % X and Mu as row vectors if they were both column vectors.
    if (d == 1) && (numel(X) > 1) && (sz(1) == n)
        X0 = X0';
        [n,d] = size(X0);
    end
    
    % Data and mean are a single row, rep them out to match covariance
    if n == 1 % already know size(Sigma,3) > 1
        n = sz(3);
        X0 = repmat(X0,n,1); % rep centered data out to match cov
    end

    % Make sure Sigma is the right size
    if sz(1) ~= sz(2)
        error(message('stats:mvggd:BadCovarianceMultiple'));
    elseif (sz(1) ~= d) || (sz(2) ~= d) % Sigma is a stack of dxd matrices
        error(message('stats:mvggd:CovSizeMismatchMultiple'));
    elseif sz(3) ~= n
        error(message('stats:mvggd:CovSizeMismatchPages'));
    else
        if sigmaIsDiag
            if any(any(Sigma<=0))
                error(message('stats:mvggd:BadDiagSigma'));
            end
            R = sqrt(Sigma);
            xRinv = X0./R;
        else
            % Create array of standardized data, and vector of log(sqrt(det(Sigma)))
            xRinv = zeros(n,d,'like',internal.stats.dominantType(X0,Sigma));
            for i = 1:n
                % Make sure Sigma is a valid covariance matrix
                [R,err] = cholcov(Sigma(:,:,i),0);
                if err ~= 0
                    error(message('stats:mvggd:BadMatrixSigmaMultiple'));
                end
                xRinv(i,:) = X0(i,:) / R;
            end
        end
    end
   
elseif ndims(Sigma) > 3
    error(message('stats:mvggd:BadCovariance'));
end

% Assume unitary Beta
if nargin < 4 || isempty(Beta)
    Beta = 1; 
    B = gamma((d+2)/(2*Beta))/(d*gamma(d/(2*Beta)));
    
    if ndims(Sigma) == 3
        A = (B/pi)^(d/2)*(gamma(d/2)*Beta)/(gamma(d/(2*Beta)).*arrayfun(@(k) det(Sigma(:, :, k)), 1:size(Sigma, 3)).^(-1/2));
    else
        A = (B/pi)^(d/2)*(gamma(d/2)*Beta)/(gamma(d/(2*Beta))*det(Sigma)^(-1/2));
    end
% Get vector Beta
elseif ismatrix(Beta)
    [n3,~] = size(Beta);
    if (n3 ~= n) && (n3 ~= 1) % has to have same number of Betas as X
        error(message('stats:mvggd:rowSizeMismatch'));
    else
        if any(Beta<=0) || any(isnan(Beta)) || any(isinf(Beta))
            error(message('stats:mvggd:BadBeta'));
        else
            B = gamma((d+2)./(2.*Beta))./(d.*gamma(d./(2.*Beta)));
            if ndims(Sigma) == 3
                A = (B./pi).^(d/2)*(gamma(d/2).*Beta)./(gamma(d./(2.*Beta)).*arrayfun(@(k) det(Sigma(:, :, k)), 1:size(Sigma, 3)).^(-1/2));
            else
                A = (B./pi).^(d/2).*(gamma(d/2).*Beta)./(gamma(d./(2.*Beta)).*det(Sigma)^(-1/2));
            end
        end
    end
else
    error(message('stats:mvggd:BadBeta'));
end

yBeforeNaNCheck = A.*exp(-(B.*sum((xRinv.^2), 2)).^ Beta);

% Special case correction: X or Mu contain Inf, but not both
% g2835301 
% NaNs in the solution, y, are corrected to 0 as needed.
% Note: In the case Sigma has an Inf, we leave the solution unchanged
% Note: We do not handle the case where Inf exists in Mu and X
% Note: If X, Mu, or Sigma have a NaN, we leave the solution unchanged
nanInY = isnan(yBeforeNaNCheck); % Location of NaN values in solution
nanInSolution = any(nanInY); % Check if there is a NaN in the solution
if nanInSolution % If yes, we need to check if the NaN should be 0
    infInXRow = any(isinf(X),2); % Find rows with Inf in X
    nanInXRow = any(isnan(X),2); % Find rows with NaN in X
    switch nargin % Based on number of function inputs
        case 1 % only X
            yBeforeNaNCheck(infInXRow & ~nanInXRow) = 0; % Replace incorrect NaNs with 0
            y = yBeforeNaNCheck;
        case 2 % X and Mu
            y = adjustSolutionGivenMu(yBeforeNaNCheck, Mu, infInXRow, nanInXRow);
        case 3 || 4 % X, Mu, and Sigma or X, Mu, Sigma, and Beta
            if ~any(isinf(Sigma), 'all') && ~any(isnan(Sigma), 'all') % If Sigma has an Inf or NaN, we will leave solution unchanged
               y = adjustSolutionGivenMu(yBeforeNaNCheck, Mu, infInXRow, nanInXRow);
            else
               y = yBeforeNaNCheck;
            end
    end
else
    y = yBeforeNaNCheck;
end
end

% Checks Mu for Inf and NaN, also considers X and updates solution if needed
function y = adjustSolutionGivenMu(y, Mu, infInXRow, nanInXRow)
    infInMuRow = any(isinf(Mu),2); % Find rows with Inf in Mu
    nanInMuRow = any(isnan(Mu),2); % Find rows with NaN in Mu
    y(xor(infInXRow, infInMuRow) & ~nanInXRow & ~nanInMuRow) = 0; % Replace incorrect NaNs with 0
end

