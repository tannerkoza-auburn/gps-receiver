function out = sdCarr3D(obj, carr_user, carr_base, svPos, basePos, N_ab)
%{
    Single-difference DGPS solution with carrier measurements. Assumes the
    integer ambiguities have already been found (Lambda method function
    provided in gps-lab3 workspace)

    Inputs:
    - carr_user: user carrier measurements (m)
    - carr_base: base carrier measurements (m)
    - svPos: satellite ECEF positions (m)
    - basePos: known (or estimated) ECEF base position (m)
    - N_ab: integer ambiguity difference between the two (user_int - base_int) 

    Output
    - out.pos: ECEF position solution of user (m)
    - out.clock_bias: clock bias btween user and base (user - base) (m)
    - out.DOP: dilution of precision 
    - out.P: solution covariance matrix

    
%}

%% Initialization

    % Handle Input Dimensions
    [carr_user, carr_base, svPos] = obj.dimHandle(carr_user, carr_base, svPos);

    % Define Number of Measurements
    numMeas = length(carr_user);

    % Initialize Shared Variables
    uhat_x = 0;
    uhat_y = 0;
    uhat_z = 0;
    y = 0;
    G = 0;

    
    
%% Estimation 

    unitVecs; % construct unit vecs
    measVecs; % construct measurment vector
    geomMatrix; % construct geometry matrix
    
    est = ( G' * G )^-1 * G' * y;
    
    % User-Base Relative Position Vector
    rpv = -est(1:3); % estimate comes out as vector from a to b... reverse it
    
    % global position for user
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
    
        y = carr_user - carr_base - N_ab;
    
    end
    
    function geomMatrix
    
    % Geometry Matrix Population
     G = [uhat_x uhat_y uhat_z ones(numMeas,1)];

    end




end

