function out = sdp3D(obj, psr_user, psr_base, svPos, basePos)
% DESCRIPTION: sdp3D produces a 3D state solution from 3D receiver and base
% station data.
%
% NOTE: Transpose svPos to where corresponding satellite position vectors
% (x,y) are in columns. Also, if basePos is unknown, run p3D using the base
% pseudoranges before running this function.
% 
% PARAMS:
%   - psr_user: user pseudoranges (m)
%   - psr_base: base station pseudoranges (m)
%   - svPos: satellite ECEF positions (m)
%   - basePos: known ECEF base position (m)
%
% OUTPUT:
%   - out.pos: ECEF position solution (m)
%   - out.clock_bias: clock bias solution (m)
%   - out.DOP: dillution of precision
%   - out.P: solution covariance matrix
%
% AUTHOR: Tanner Koza

%% Initialization

    % Handle Input Dimensions
    [psr_user, psr_base, svPos] = obj.dimHandle(psr_user, psr_base, svPos);

    % Define Number of Measurements
    numMeas = length(psr_user);

    % Initialize Shared Variables
    uhat_x = 0;
    uhat_y = 0;
    uhat_z = 0;
    y = 0;
    G = 0;

%% Estimation

    unitVecs

    measVec

    geomMatrix

    est = ( G' * G )^-1 * G' * y;

    rpv = -est(1:3);

    % User-Base Relative Position Vector
    pos = basePos + rpv;

    % Covariance & DOP
    DOP = ( G' * G )^-1;
    P = obj.rcvrSigma * DOP;

%% Solution Structure Population

    % Populate Structure
    out.pos = pos;
    out.rpv = rpv;
    out.clock_bias = est(4);
    out.DOP = DOP;
    out.P = P;

%% Nested Functions

    function unitVecs
    
        % Initialization
        uhat_x = zeros(numMeas,1);
        uhat_y = zeros(numMeas,1);
        uhat_z = zeros(numMeas,1);
    
        % Calculate Satellite Unit Vectors
        for i = 1:numMeas
    
            r = sqrt( ( svPos(1,i) - basePos(1) )^2 ...
                + ( svPos(2,i) - basePos(2) )^2 ...
                + ( svPos(3,i) - basePos(3) )^2 );  
    
            uhat_x(i) = ( svPos(1,i) - basePos(1) )/ r;
    
            uhat_y(i) = ( svPos(2,i) - basePos(2) )/ r;

            uhat_z(i) = ( svPos(3,i) - basePos(3) )/ r;
    
        end
    
    end
    
    function measVec
    
        y = psr_user - psr_base;
    
    end
    
    function geomMatrix
    
    % Geometry Matrix Population
     G = [uhat_x uhat_y uhat_z ones(numMeas,1)];

    end

end