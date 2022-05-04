function out = p3D(obj, psr, svPos, svClockCorr)
% DESCRIPTION: pv3D produces a GNSS state solution from GNSS receiver data.
%
% PARAMS:
%   - psr: pseudoranges (m)
%   - svPos: satellite ECEF positions (m)
%   - svClockCorr: satellite clock corrections (s)

%
% OUTPUT:
%   - out.pos: ECEF position solution (m)
%   - out.vel: ECEF velocity solution (m/s)
%   - out.clock_bias: clock bias solution (m)
%   - out.clock_drift: clock drift solution (m/s)
%   - out.DOP: dillution of precision
%   - out.P: solution covariance matrix
%
% AUTHOR: Tanner Koza

%% Initialization

% Handle Input Dimensions
[psr, svPos] = obj.dimHandle(psr, svPos);

% Define Number of Measurements
numMeas = length(psr);

% Initialize Shared Variables
uhat_x = 0;
uhat_y = 0;
uhat_z = 0;
y = 0;
G = 0;

% Initialize State Estimate Vector
est = [obj.initPos; obj.initClockBias];

% Initialize Least Squares Iteration Count
itr = 0;

% SV Clock Correction Unit Conversion (s to m)
C = physconst('LightSpeed');
svClockCorr = svClockCorr * C;

%% Estimation

% Least Squares & Newton-Raphson
while true

    unitVecs

    measVec

    geomMatrix

    dest = ( G' * G )^-1 * G' * y;

    est = est + dest;

    itr = itr + 1;

    if  norm( dest(1:4) ) <= obj.lsconv

        break

    end

end

% Covariance & DOP
DOP = ( G' * G )^-1;
P = obj.rcvrSigma * DOP;

%% Solution Structure Population

% Populate Structure
out.pos = est(1:3);
out.clock_bias = est(4);
out.DOP = DOP;
out.P = P;

% Reinitialize Class Variables
% NOTE: Reinitializing initVel and initCloclDrift are unecessary
% because their final values are not invloved in Newton-Raphson.
obj.initPos = est(1:3);
obj.initClockBias = est(4);

%% Nested Functions

    function unitVecs

        % Initialization
        uhat_x = zeros(numMeas,1);
        uhat_y = zeros(numMeas,1);
        uhat_z = zeros(numMeas,1);

        % Calculate Satellite Unit Vectors
        for i = 1:numMeas

            r = sqrt( ( svPos(1,i) - est(1) )^2 ...
                + ( svPos(2,i) - est(2) )^2 ...
                + ( svPos(3,i) - est(3) )^2);

            uhat_x(i) = ( svPos(1,i) - est(1) )/ r;

            uhat_y(i) = ( svPos(2,i) - est(2) )/ r;

            uhat_z(i) = ( svPos(3,i) - est(3) )/ r;

        end

    end

    function measVec

        % Initialization
        psrhat = zeros(numMeas,1);

        % Measurement Vector Population
        for i = 1:numMeas

            if svClockCorr(i) == 0
                psrhat(i) = sqrt( ( svPos(1,i) - est(1) )^2 ...
                    + ( svPos(2,i) - est(2) )^2 ...
                    + ( svPos(3,i) - est(3) )^2) - svClockCorr(i);
            else
                psrhat(i) = sqrt( ( svPos(1,i) - est(1) )^2 ...
                    + ( svPos(2,i) - est(2) )^2 ...
                    + ( svPos(3,i) - est(3) )^2) + est(4) - svClockCorr(i);
            end

        end

        y = psr - psrhat;

    end

    function geomMatrix

        % Geometry Matrix Population
        G = [-uhat_x -uhat_y -uhat_z ones(numMeas,1)];

    end

end