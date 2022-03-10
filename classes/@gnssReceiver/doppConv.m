function out = doppConv(obj, dopp, carrFreq)
% DESCRIPTION: doppConv converts GNSS doppler measurements (Hz) into speed
% measurements (m/s) based on a particular carrier frequency.
% PARAMS:
%   - dopp: GNSS doppler measurements (Hz)
%   - carrFreq: GPS carrier L band
% OUTPUT:
%   - out: GNSS doppler measurements (m/s)
%
% AUTHOR: Tanner Koza


%% Initialization

    % Define Speed of Light Constant
    C = physconst('LightSpeed');

%% Doppler Conversion

    if(carrFreq == 1)
        
        out = dopp * -( C/obj.L1 );

    elseif(carrFreq == 2)

        out = dopp * -( C/obj.L2 );

    elseif(carrFreq == 3)

        out = dopp * -( C/obj.L3 );

    elseif(carrFreq == 5)

        out = dopp * -( C/obj.L5 );

    end

end