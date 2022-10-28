function out = ddp3D(obj, psr_user, psr_base, svPos, basePos)
% DESCRIPTION: ddp2D produces a 2D state solution from 2D receiver and base
% station data.
%
% NOTE: Transpose svPos to where corresponding satellite position vectors
% (x,y) are in columns. Also, if basePos is unknown, run p2D using the base
% pseudoranges before running this function.
% 
% PARAMS:
%   - psr_user: user pseudoranges (m)
%   - psr_base: base station pseudoranges (m)
%   - svPos: satellite ECEF positions (m)
%   - basePos: known ECEF base position (m)
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
    uhat_z = 0;l
    y = zeros(numMeas-1,1);
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
    
        delta_psr = psr_user - psr_base;

        for i = 2:numMeas

            y(i-1) = delta_psr(1) - delta_psr(i);

        end

    
    end
    
    function geomMatrix
    
        delta_uhat_x = zeros(numMeas-1,1);
        delta_uhat_y = zeros(numMeas-1,1);
        delta_uhat_z = zeros(numMeas-1,1);

    for i = 2:numMeas

        delta_uhat_x(i-1) = uhat_x(1) - uhat_x(i);
        delta_uhat_y(i-1) = uhat_y(1) - uhat_y(i);
        delta_uhat_z(i-1) = uhat_z(1) - uhat_z(i);

    end
    % Geometry Matrix Population
     G = [delta_uhat_x delta_uhat_y delta_uhat_z];

    end

end