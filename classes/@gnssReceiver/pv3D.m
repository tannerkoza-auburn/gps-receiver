function out = pv3D(obj, psr, dopp, svPos, svVel, svClockCorr, carrFreq)
% DESCRIPTION: pv3D produces a GNSS state solution from GNSS receiver data.
%
% PARAMS:
%   - psr: pseudoranges (m)
%   - dopp: doppler measurements (Hz)
%   - svPos: satellite ECEF positions (m)
%   - svVel: satellite ECEF velocities (m/s)
%   - svClockCorr: satellite clock corrections (s)
%   - carrFreq: GPS carrier L band
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
    [psr, dopp, svPos, svVel] = obj.dimHandle(psr, dopp, svPos, svVel);

    % Define Number of Measurements
    numMeas = length(psr);

    % Initialize Shared Variables
    uvs = 0;
    uhat_x = 0;
    uhat_y = 0;
    uhat_z = 0;
    y = 0;
    G = 0;

    % Initialize State Estimate Vector
    est = [obj.initPos; obj.initClockBias; obj.initVel; obj.initClockDrift];
    
    % Initialize Least Squares Iteration Count
    itr = 0;
    
    % Doppler Measurement Unit Conversion (Hz to m/s)
    psrdot = obj.doppConv(dopp, carrFreq);
    
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
    out.vel = dest(5:7);
    out.clock_bias = est(4);
    out.clock_drfit = dest(8);
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
    
        uvs = [uhat_x uhat_y uhat_z];
    
    end
    
    function measVec
    
        % Initialization 
        psrhat = zeros(numMeas,1);
        psrdothat = zeros(numMeas,1);
    
        % Measurement Vector Population
        for i = 1:numMeas
    
            psrhat(i) = sqrt( ( svPos(1,i) - est(1) )^2 ...
                + ( svPos(2,i) - est(2) )^2 ...
                + ( svPos(3,i) - est(3) )^2) + est(4) - svClockCorr(i);
    
            psrdothat(i) = uvs(i,1) * svVel(1,i) + ...
                           uvs(i,2) * svVel(2,i) + uvs(i,3) * svVel(3,i);
    
        end
    
        y = [psr; psrdot] - [psrhat; psrdothat];
    
    end
    
    function geomMatrix
    
    % Geometry Matrix Population
     G = [-uhat_x -uhat_y -uhat_z ones(numMeas,1) zeros(numMeas,4);
        zeros(numMeas,4) -uhat_x -uhat_y -uhat_z ones(numMeas,1)];

    end

end