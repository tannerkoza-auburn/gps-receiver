classdef gnssReceiver < handle

    properties (Access = public)

        % Receiver Properties
        rcvrSigma

        % State Estimate Initialization
        initPos = [0; 0; 0]; % ECEF (m)
        initVel = [0; 0; 0]; % ECEF (m/s)
        initClockBias = 0; % (m)
        initClockDrift = 0; % (m/s)

    end

    properties (Access = protected, Constant)

        % Carrier Frequency Definitions (Hz)
        L1 = 1575.42e6;
        L2 = 1227.60e6;
        L3 = 1381.05e6;
        L5 = 1176.45e6;

        % Least Squares Convergence Criterion
        lsconv = 1e-6;

    end

    methods (Access = public)

        % Class Constructor
        function obj = gnssReceiver(rcvrSigma)

            if nargin == 0
                
                % Default Receiver Accuracy (m)
                obj.rcvrSigma = 0.5;

            else
                
                obj.rcvrSigma = rcvrSigma;

            end
            
        end
        
        % 3D Position Estimate
        p3D(obj)
        
        % 3D Position & Velocity Estimate
        est = pv3D(obj, psr, dopp, svPos, svVel, svClockCorr, carrFreq)
        
        % 2D Position Estimate
        est = p2D(obj, psr, svPos)

        % 2D Position Estimate with Perfect Clock
        est = p2DPC(obj, psr, svPos)

        % 2D Signle Difference Position Estimate
        est = sdp2D(obj, psr_user, psr_base, svPos, basePos)
                
    end

    methods (Access = protected, Hidden, Sealed)

        % Doppler Conversion
        dopp = doppConv(obj, dopp, carrFreq)

    end

    methods (Access = protected, Hidden, Sealed, Static)
        
        % Input Data Dimension Handling
        varargout = dimHandle(varargin)

    end

end


