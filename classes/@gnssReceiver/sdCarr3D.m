function out = sdCarr3D(obj, psr_user, carrL1_user, carrL2_user, psr_base, carrL1_base, carrL2_base, svPos, basePos)
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

    
    
    there exists an integer ambiguity btwn base and user for each satellite and 
    each frequency... eg: 6 SVs = 12 ints to solve lambdaMat is numSats*3 X numSats*2

%}

%% Initialization

    % Handle Input Dimensions
    [psr_user, carrL1_user, carrL2_user, psr_base, carrL1_base, carrL2_base, svPos] = ...
        obj.dimHandle(psr_user, carrL1_user, carrL2_user, psr_base, carrL1_base, carrL2_base, svPos);
    
    % Define Number of Measurements (3 * SVs for float solution)
    numMeas = length(carrL1_user);
    
    % Initialize Shared Variables
    uhat_x = 0;
    uhat_y = 0;
    uhat_z = 0;
    yN = 0;
    y = [];
    G = [];
    G_N = [];
    lambdaMat = [];
    
    
%% First construct the float-solution for the integer ambiguities
%{ 
going to use left-null to calculate Ly = L*G*r + L*lambda*N with least squares
measurements ordered as [rho_1, phi_L1_1, phi_L2_1, rho_2, phi_L1_2, phi_L2_2...
3x SVs measurements for the float solution
%}    
    unitVecs; % construct unit vecs
    geomMatrixN; % construct geometry matrix (G) and the lambda matrix - same form for the psr and carr
    
    % extend lambda matrix and construct y for number of sats in view:    
    for j = 1:numMeas        
        
        delRho = psr_user(j) - psr_base(j);
        delPhi_L1 = carrL1_user(j) - carrL1_base(j);
        delPhi_L2 = carrL2_user(j) - carrL2_base(j);
        
        if j == 1
            yN = [delRho; delPhi_L1; delPhi_L2];
            R = [1^2    0   0; 
                0   .01^2     0;
                0   0   .01^2]; % noise on psr, phi, phi
        else
            yN = [yN; delRho; delPhi_L1; delPhi_L2];
            R = blkdiag(R, [1^2    0   0; 
                        0   .01^2     0;
                        0   0   .01^2]);
        end         
    end 
    
    
    %--- float solution of integer ambiguities
    L = (null(G_N.'))'; % calculate left Null: L*G = 0
    
    % remap matrices with left-null for least-squares form
    H = L*lambdaMat;
    yN = L*yN;
    R = L*R*L.'; % not in notes... pulled this out of thin air based on the sizing constraints for matrix algebra
    
    % float-level integer estimates
    N_est = (H.' * H)^(-1)*H.'*yN;
    
    % float-level integer uncertainty
    P_N = inv(H.' * inv(R) * H);
    
    
%% Lambda Corrections    
    
    
    
%% Estimation 

    measVec; % construct measurement vector with new N estimates (L1 carriers only)
    geomMatrix; % create L1 carrier geometry matrix
    
    est = ( G' * G)^-1 * G' * y;
    
    % User-Base Relative Position Vector
    rpv = -est(1:3); % estimate comes out as vector from a to b... reverse it
    
    % global position for user
    pos = basePos + rpv;
    
    % Covariance & DOP
    DOP = ( G' * G )^-1;
    P = 0.01^2 * DOP; % assumed carrier sigma

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
    
        for i = 1:numMeas
            N_L1 = N_est(1 + 2*(i - 1)); % odd indices are L1 integers
            y(i,1) = carrL1_user(j) - carrL1_base(j) - N_L1;
        end 
    
    end
    
    function geomMatrix
        % construct Phi_L1 measurement matrix... not using psr or L2    
        for i = 1:numMeas            
            G(end+1,:) = [uhat_x(i) uhat_y(i) uhat_z(i) 1];
        end 
    end 

    function geomMatrixN
    
    % Geometry and lambda Matrix Population
    C = physconst('LightSpeed');
    lamb1 = C / 1575.42*10^6;
    lamb2 = C / 1227.6*10^6;
    
    for i = 1:numMeas
        
        % 3x for each sat... delRho, delPhiL1, delPhiL2
        G_N(end+1:end+3, :) = [uhat_x(i) uhat_y(i) uhat_z(i) 1;
                             uhat_x(i) uhat_y(i) uhat_z(i) 1;
                             uhat_x(i) uhat_y(i) uhat_z(i) 1];       
                
        if i == 1
            lambdaMat =  [0 0;
                        lamb1 0; 
                        0   lamb2];
        else 
            lambdaMat = blkdiag(lambdaMat, [0 0;
                                            lamb1 0; 
                                            0   lamb2]);
        end 
    end 

    end




end

