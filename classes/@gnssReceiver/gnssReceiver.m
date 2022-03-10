classdef gnssReceiver < handle
% DESCRIPTION: gnssReceiver is a class that contains methods to produce
% different state solutions depending on the input data to each method.
% There are 2D and 3D solutions that assume receiver clock bias and a
% perfect clock. DGPS solutions are also available.
%
% METHODS: 
%   - p3D: 3 dimensional position solution
%   - pv3D: 3 dimensional position & velocity solution
%   - p2D: 2 dimensional position solution
%   - p2DPC: 2 dimensional position solution assuming perfect clock
%   - sdp2D: single difference 2 dimensional position solution
%   - ddp2D: double difference 2 dimensional position solution
%
% AUTHOR: Tanner Koza

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

        % 2D Single Difference Position Estimate
        est = sdp2D(obj, psr_user, psr_base, svPos, basePos)

        % 2D Double Difference Position Estimate
        est = ddp2D(obj, psr_user, psr_base, svPos, basePos)
                
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


