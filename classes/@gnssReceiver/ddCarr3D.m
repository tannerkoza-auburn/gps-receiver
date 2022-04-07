function out = sdCarr3D(obj, psr_user, carrL1_user, carrL2_user, psr_base, carrL1_base, carrL2_base, svPos, basePos)
%{
    Single-difference DGPS solution with carrier measurements. Solves for
    the integer ambiguities with a widelane approach using L1 and L2 toghether.

    Inputs:
    - psr_user: user pseudoranges (L1)
    - carr_user: user carrier measurements (L1 and L2) (m)
    - psr_base: base pseudoranges (L1)
    - carr_base: base carrier measurements (L1 and L2) (m)
    - svPos: satellite ECEF positions (m)
    - basePos: known (or estimated) ECEF base position (m)

    Output
    - out.pos: ECEF position solution of user (m)
    - out.N_lambda: resulting lambda integer estimates
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
    for j = 2:numMeas        
        
        dblDelRho = (psr_user(1) - psr_base(1)) - (psr_user(j) - psr_base(j));
        dblDelPhi_L1 = (carrL1_user(1) - carrL1_base(1)) - (carrL1_user(j) - carrL1_base(j));
        dblDelPhi_L2 = (carrL2_user(1) - carrL2_base(1)) - (carrL2_user(j) - carrL2_base(j));
        
        if j == 2
            yN = [dblDelRho; dblDelPhi_L1; dblDelPhi_L2];
            R = [sqrt(1^2 + 1^2)    0   0; 
                0   sqrt(.01^2 + .01^2)     0;
                0   0   sqrt(.01^2 + 0.1^2)]; % noise on psr, phi, phi
        else
            yN = [yN; dblDelRho; dblDelPhi_L1; dblDelPhi_L2];
            R = blkdiag(R, [sqrt(1^2 + 1^2)    0   0; 
                0   sqrt(.01^2 + .01^2)     0;
                0   0   sqrt(.01^2 + 0.1^2)]); % noise on psr, phi, phi
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
    P_N = (P_N + P_N.')/2; % force symmetry! 
        
%% Lambda Corrections    

    [N_est, sqnorm] = LAMBDA(N_est, P_N, 1);
    N_est = N_est(:,1);
    out.N_lambda = N_est; % save estimate of integers to ouput struct
    % ^ above reports mutliple options for integers along with estimates of
    % accuracy for each... 
    
    
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
    out.DOP = DOP;
    out.P = P;

%% Nested Functions

     function unitVecs % calculate delta unit vectors
    
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
            C = physconst('LightSpeed');
            lamb1 = C / (1575.42e6);
            %lamb2 = C / (1227.6e6);
           % N_L1 = N_est(2*i - 1); % odd indices are L1 integers
            %N_L2 = N_est(2*i); % even indices are L2 integers
            
            delCarrL1(i) = carrL1_user(i) - carrL1_base(i);
            
            % y(i,1) = carrL2_user(i) - carrL2_base(i) - lamb2*N_L2; % L2 option
        end 
    
        for i = 2:numMeas
            N_L1 = N_est(2*i - 3); % odd indices are L1 integers
            y(i-1,1) = delCarrL1(1) - delCarrL1(i) - lamb1*N_L1; % dbl difference
        end
            
        
    end
    
    function geomMatrix
        % construct Phi_L1 measurement matrix... not using psr or L2    
        for i = 2:numMeas            
            G(end+1,:) = G_N(3*(i-1),:); %[uhat_x(i) uhat_y(i) uhat_z(i) 1];
        end 
    end 

    function geomMatrixN
    
    % Geometry and lambda Matrix Population
    C = physconst('LightSpeed');
    lamb1 = C / (1575.42e6);
    lamb2 = C / (1227.6e6);
    

    for i = 2:numMeas
        % delta unit vectors between first sat and others
        delta_uhat_x = uhat_x(1) - uhat_x(i);
        delta_uhat_y = uhat_y(1) - uhat_y(i);
        delta_uhat_z = uhat_z(1) - uhat_z(i);
        
        if i == 2
            lambdaMat =  [0 0;
                        lamb1 0; 
                        0   lamb2];
        else 
            lambdaMat = blkdiag(lambdaMat, [0 0;
                                            lamb1 0; 
                                            0   lamb2]);
        end 
        
        G_N(end+1:end+3, :) = [delta_uhat_x delta_uhat_y delta_uhat_z;
                             delta_uhat_x delta_uhat_y delta_uhat_z;
                             delta_uhat_x delta_uhat_y delta_uhat_z];    
        
    end 

    end




end

