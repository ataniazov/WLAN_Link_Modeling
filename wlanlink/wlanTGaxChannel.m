classdef (StrictDefaults)wlanTGaxChannel < wlan.internal.ChannelBase
%wlanTGaxChannel Filter input signal through a TGax multipath fading channel
%   tgax = wlanTGaxChannel creates a System object, tgax, for the TGax
%   indoor MIMO channel model as specified by the IEEE 802.11 Wireless LAN
%   Working group [1,2,3,4], which follows the MIMO modeling approach
%   presented in [2]. This object filters an input signal through the
%   multipath TGax channel to obtain the channel impaired signal.
%
%   tgax = wlanTGaxChannel(Name,Value) creates a TGax channel object, tgax,
%   with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   Step method syntax:
%
%   Y = step(tgax,X) filters input signal X through a TGax fading channel
%   and returns the result in Y. The input X must be a double precision,
%   real or complex matrix of size Ns-by-Nt where Ns is the number of
%   samples and Nt is the number of transmit antennas. Nt must be equal to
%   the NumTransmitAntennas property value of tgax. Y is the output signal
%   of size Ns-by-Nr, where Nr is the number of receive antennas that takes
%   the same value as the NumReceiveAntennas property value of tgax. Y is
%   of double precision data type with complex values.
% 
%   [Y,PATHGAINS] = step(tgax,X) returns the TGax channel path gains of the
%   underlying fading process in PATHGAINS. This syntax applies when you
%   set the PathGainsOutputPort property of tgax to true. PATHGAINS is of
%   size Ns-by-Np-by-Nt-by-Nr, where Np is the number of resolvable paths,
%   that is, the number of paths defined for the case specified by the
%   DelayProfile property. PATHGAINS is of double precision data type with
%   complex values.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(obj,x) and y = obj(x) are
%   equivalent.
% 
%   wlanTGaxChannel methods:
%
%   step     - Filter the input signal through a MIMO fading channel
%   release  - Allow change of property values and input characteristics
%   clone    - Create TGax channel object with same property values
%   isLocked - Locked status (logical)
%   <a href="matlab:help matlab.System/reset   ">reset</a>    - Reset state of filters and random stream if the
%              RandomStream property is set to 'mt19937ar with seed'
%   <a href="matlab:help matlab.System/info   ">info</a>     - Return characteristic information about the TGax channel 
%
%   wlanTGaxChannel properties:
%
%   SampleRate              - Input signal sample rate (Hz)
%   DelayProfile            - Delay profile models for WLAN
%   ChannelBandwidth        - Channel bandwidth of the input waveform
%   CarrierFrequency        - Carrier frequency (Hz)
%   EnvironmentalSpeed      - Speed of the scatterers (km/h)
%   TransmitReceiveDistance - Distance between transmit and receive (m)
%   NormalizePathGains      - Normalize path gains (logical)
%   UserIndex               - UserIndex for single and multiuser scenario
%   TransmissionDirection   - Transmission direction (Uplink/Downlink)
%   NumTransmitAntennas     - Number of transmit antennas
%   TransmitAntennaSpacing  - Transmit antenna spacing in wavelength
%   NumReceiveAntennas      - Number of receive antennas
%   ReceiveAntennaSpacing   - Receive antenna spacing in wavelength
%   LargeScaleFadingEffect  - Inclusion of large scale fading effect
%   NumPenetratedFloors     - Number of building floors between transmitter and receiver
%   NumPenetratedWalls      - Number of walls between transmitter and receiver
%   WallPenetrationLoss     - Penetration loss of a single wall (dB)
%   FluorescentEffect       - Enable fluorescent effect in channel modeling (logical)
%   PowerLineFrequency      - Power line frequency (Hz)
%   NormalizeChannelOutputs - Normalize channel outputs (logical)
%   RandomStream            - Source of random number stream
%   Seed                    - Initial seed of mt19937ar random number stream
%   PathGainsOutputPort     - Enable path gain output (logical)
%
%   % Example: 
%   %   Get a channel impulse response by filtering an impulse through a
%   %   TGax channel.
%
%   input = zeros(100,1);
%   input(10) = 1; % Create an impulse
%   tgax = wlanTGaxChannel;
%   tgax.LargeScaleFadingEffect = 'PathLoss and shadowing';
%   tgax.NumPenetratedFloors = 2;
%   tgax.RandomStream = 'mt19937ar with seed';
%   tgax.Seed = 10;
%   tgax.SampleRate = 1e9;
%   figure;
%   time = (1/tgax.SampleRate)*(0:length(input)-1);
%   stem(time,abs(tgax(input)));
%   xlabel('Time (sec)');
%   ylabel('|Amplitude|');
%   title('Channel Impulse Response');
%
%   % References:
%   % [1] Jianhan. L, Ron. P, et al. TGax Channel Model. IEEE
%   % 802.11-14/0882r4, September 2014.
%   % [2] Erceg, V., L. Schumacher, P. Kyritsi, et al. TGn Channel Models.
%   % Version 4. IEEE 802.11-03/940r4, May 2004.
%   % [3] Breit, G., H. Sampath, S. Vermani, et al. TGac Channel Model
%   % Addendum. Version 12. IEEE 802.11-09/0308r12, March 2010.
%   % [4] Kermoal, J. P., L. Schumacher, K. I. Pedersen, P. E. Mogensen, and
%   % F. Frederiksen, "A Stochastic MIMO Radio Channel Model with
%   % Experimental Validation", IEEE Journal on Selected Areas in
%   % Communications, Vol. 20, No. 6, August 2002, pp. 1211-1226.
%
%   See also wlanTGnChannel, wlanTGacChannel, wlanTGahChannel.

% Copyright 2017 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

% Public properties
properties (Nontunable)
    %SampleRate Sample rate (Hz)
    %   Specify the sample rate of the input signal in Hz as a double
    %   precision, real, positive scalar. The default is 80e6 Hz.
    SampleRate = 80e6;
    %ChannelBandwidth Channel bandwidth of the input waveform
    %   Specify the channel bandwidth of the input waveform as one of
    %   'CBW20'|'CBW40'|'CBW80'|'CBW160'. The ChannelBandwidth is used to
    %   reduce the multipath spacing of the power delay profile by a factor
    %   of 2^ceil(log2(BW/40)), where BW is the channel bandwidth in MHz.
    %   The reduction factor applies for channel bandwidths greater than 40
    %   MHz. The default is 'CBW80'.
    ChannelBandwidth = 'CBW80'
    %UserIndex User index for single or multiuser scenario
    %   Specify a particular user in a multiuser scenario. A pseudorandom
    %   per-user angle-of-arrival (AoA) and angle-of-departure (AoD)
    %   rotation is applied to support multiuser scenario. The value of
    %   zero means a simulation scenario not requiring per user angle
    %   diversity and assuming the TGn defined cluster AoAs and AoDs. The
    %   default is zero.
    UserIndex = 0;
    %TransmissionDirection Transmission direction (Uplink/Downlink)
    %   Specify the transmission direction as one of 'Uplink' | 'Downlink'.
    %   The default is 'Downlink'.
    TransmissionDirection = 'Downlink';
    %NumTransmitAntennas Number of transmit antennas
    %   Specify the number of transmit antennas as a numeric, real,
    %   positive integer scalar between 1 and 8, inclusive. The default is
    %   1.
    NumTransmitAntennas = 1;
    % TransmitAntennaSpacing Transmit antenna spacing in wavelengths
    %   Spacing of the regular geometry of the antenna elements at the
    %   transmitter, in wavelengths. Only uniform linear array is
    %   supported. This property applies only when NumTransmitAntennas is
    %   greater than 1. The default is 0.5.
    TransmitAntennaSpacing = 0.5;
    %NumReceiveAntennas Number of receive antennas
    %   Specify the number of receive antennas as a numeric, real, positive
    %   integer scalar between 1 and 8, inclusive. The default is 1.
    NumReceiveAntennas = 1;
    % ReceiveAntennaSpacing Receive antenna spacing in wavelengths
    %   Spacing of the regular geometry of the antenna elements at the
    %   receiver, in wavelengths. Only uniform linear array is supported.
    %   This property applies only when NumReceiveAntennas is greater than
    %   1. The default is 0.5.
    ReceiveAntennaSpacing = 0.5;
    % CarrierFrequency Carrier frequency (Hz) 
    %   Specify the carrier frequency of the input signal in Hz. The
    %   default is 5.25e9 Hz.
    CarrierFrequency = 5.25e9;
    % EnvironmentalSpeed Speed of the scatterers (km/h)
    %   Specify the speed of the scatterers in km/h as numeric positive
    %   scalar greater than or equal to zero. The TGn and TGac channel
    %   documents specify an environment speed of 1.2 km/h and 0.089 km/h
    %   respectively. The default is 0.089 km/h.
    EnvironmentalSpeed = 0.089;
    % NumPenetratedFloors Number of building floors
    %   Specify the number of building floors between the transmitter and
    %   the receiver to account for the floor attenuation loss in the
    %   calculation of path loss in a multiple floor scenario. The default
    %   is 0, which represents a communication link between a transmitter
    %   and a receiver located on the same floor.
    NumPenetratedFloors = 0;
    % NumPenetratedWalls Number of walls
    %   Specify the number of walls between the transmitter and the
    %   receiver to account for the wall penetration loss in the
    %   calculation of path loss. The default is 0, which represents a
    %   communication link between a transmitter and a receiver without
    %   wall penetration loss.
    NumPenetratedWalls = 0;
    % WallPenetrationLoss Penetration loss of a single wall (dB)
    %   Specify the penetration loss of a single wall in decibels. This
    %   property applies only when NumPenetratedWalls is greater than 0.
    %   The default value is 5dB.
    WallPenetrationLoss = 5;
end

properties(Constant,Hidden)
    ChannelBandwidthSet = matlab.system.StringSet({'CBW20','CBW40','CBW80','CBW160'});
    TransmissionDirectionSet = matlab.system.StringSet({'Downlink','Uplink'});
    InterpolationFactor = 1/40;
end

methods
  function obj = wlanTGaxChannel(varargin) % Constructor
    setProperties(obj, nargin, varargin{:});
    obj.pLegacyGenerator = false;
  end 
  
  function set.SampleRate(obj, val)
    propName = 'SampleRate';
        validateattributes(val, {'double'}, ...
            {'real','scalar','positive','finite'}, ...
            [class(obj) '.' propName], propName);   
    obj.SampleRate = val;
  end
  
  function set.UserIndex(obj, val)
    propName = 'UserIndex';
    validateattributes(val, {'numeric'}, ...
        {'real','scalar','integer','>=',0}, ...
        [class(obj) '.' propName], propName);
    obj.UserIndex= val;
  end
  
  function set.NumTransmitAntennas(obj, val)
    propName = 'NumTransmitAntennas';
    validateattributes(val, {'numeric'}, ...
        {'real','scalar','integer','>=',1,'<=',8}, ...
        [class(obj) '.' propName], propName);
    obj.NumTransmitAntennas = val;
  end
  
  function set.TransmitAntennaSpacing(obj, val)
    propName = 'TransmitAntennaSpacing';
    validateattributes(val, {'numeric'}, ...
        {'real','scalar','>',0,'finite'}, ...
        [class(obj) '.' propName], propName);
    obj.TransmitAntennaSpacing = val;
  end
  
  function set.NumReceiveAntennas(obj, val)
    propName = 'NumReceiveAntennas';
    validateattributes(val, {'numeric'}, ...
        {'real','scalar','integer','>=',1,'<=',8}, ...
        [class(obj) '.' propName], propName);
    obj.NumReceiveAntennas = val;
  end

  function set.ReceiveAntennaSpacing(obj, val)
    propName = 'ReceiveAntennaSpacing';
    validateattributes(val, {'numeric'}, ...
        {'real','scalar','>',0,'finite'}, ...
        [class(obj) '.' propName], propName);
    obj.ReceiveAntennaSpacing = val;
  end
  
  function set.CarrierFrequency(obj, val)
     propName = 'CarrierFrequency';
        validateattributes(val, {'numeric'}, ...
            {'real','scalar','>',0,'finite'}, ...
            [class(obj) '.' propName], propName); 
      obj.CarrierFrequency = val;
  end
  
  function set.EnvironmentalSpeed(obj, val)
    propName = 'EnvironmentalSpeed';
    validateattributes(val, {'numeric'}, ...
        {'real','scalar','>=',0,'finite'}, ...
        [class(obj) '.' propName], propName);
    obj.EnvironmentalSpeed = val;
  end
   
  function set.NumPenetratedFloors(obj, val)
    propName = 'NumPeneratedFloors';
    validateattributes(val, {'numeric'}, ...
        {'real','scalar','integer','>=',0,'finite'}, ...
        [class(obj) '.' propName], propName);
    obj.NumPenetratedFloors = val;
  end
  
  function set.NumPenetratedWalls(obj, val)
    propName = 'NumPeneratedWalls';
    validateattributes(val, {'numeric'}, ...
        {'real','scalar','integer','>=',0,'finite'}, ...
        [class(obj) '.' propName], propName);
    obj.NumPenetratedWalls = val;
  end
  
  function set.WallPenetrationLoss(obj, val)
    propName = 'WallPenetrationLoss';
    validateattributes(val, {'numeric'}, ...
        {'real','scalar','>=',0,'finite'}, ...
        [class(obj) '.' propName], propName);
    obj.WallPenetrationLoss = val;
  end
end

methods(Access = protected)
  function flag = isInactivePropertyImpl(obj, prop)
    % Use the if-else format for codegen
    if strcmp(prop, 'TransmitAntennaSpacing')
        flag = obj.NumTransmitAntennas == 1;
    elseif strcmp(prop, 'ReceiveAntennaSpacing')
        flag = obj.NumReceiveAntennas == 1;
    elseif strcmp(prop, 'NumPenetratedFloors') || strcmp(prop,'NumPenetratedWalls') || strcmp(prop, 'WallPenetrationLoss')
        flag = strcmp(obj.LargeScaleFadingEffect,'None') || strcmp(obj.LargeScaleFadingEffect,'Shadowing');
    else
        flag = isInactivePropertyImpl@wlan.internal.ChannelBase(obj, prop);
    end
  end

  function s = infoImpl(obj)
    % info Returns characteristic information about TGax channel
    %    S = info(OBJ) returns a structure containing characteristic
    %    information, S, about the TGax fading channel. A description of
    %    the fields and their values are as follows:
    %     
    %    ChannelFilterDelay        - Channel filter delay, measured in
    %                                samples
    %    ChannelFilterCoefficients - Channel fractional delay filter
    %                                coefficients
    %    PathDelays                - Multipath delay of the discrete paths
    %                                in seconds for TGax defined delay
    %                                profile.
    %    AveragePathGains          - Average gains of the discrete paths
    %                                in dB for TGax defined delay profile.
    %    Pathloss                  - Path loss in dBs. 
    
    if ~coder.target('MATLAB') || ~isLocked(obj)
        getInfoParameters(obj);
    end
    
    s.ChannelFilterDelay        = obj.pChannelFilterDelay;
    s.ChannelFilterCoefficients = getExpandedInterpMatrix(obj); 
    s.PathDelays                = obj.pPathDelays; 
    s.AveragePathGains          = obj.pPathPowerdBs;
    s.Pathloss                  = 0;
     
    if strcmp(obj.LargeScaleFadingEffect,'Pathloss') || ...
            strcmp(obj.LargeScaleFadingEffect,'Pathloss and shadowing')
        s.Pathloss = obj.pPathloss;
    end
  end
  
end

methods(Static, Access = protected)
  function groups = getPropertyGroupsImpl
    multipath = matlab.system.display.Section( ...
        'PropertyList',{'SampleRate','DelayProfile','ChannelBandwidth', ...
        'CarrierFrequency','EnvironmentalSpeed','TransmitReceiveDistance',  ...
        'NormalizePathGains'});
    
    antenna = matlab.system.display.Section(...
        'PropertyList', {'UserIndex','TransmissionDirection', ...
        'NumTransmitAntennas','TransmitAntennaSpacing', ...
        'NumReceiveAntennas','ReceiveAntennaSpacing'}); 
    
    pathloss = matlab.system.display.Section(...
        'PropertyList', {'LargeScaleFadingEffect','NumPenetratedFloors', ...
        'NumPenetratedWalls','WallPenetrationLoss', ...
        'FluorescentEffect','PowerLineFrequency'}); 
    
    pRandStream = matlab.system.display.internal.Property(...
        'RandomStream', ...
        'IsGraphical',false, ...
        'UseClassDefault',false, ...
        'Default','mt19937ar with seed');
    
    randomization = matlab.system.display.Section(...
        'PropertyList',{pRandStream,'Seed','PathGainsOutputPort'});
    
    normalization = matlab.system.display.Section(...
        'PropertyList',{'NormalizeChannelOutputs'});
    
    mainGroup = matlab.system.display.SectionGroup(...
        'TitleSource', 'Auto', ...
        'Sections', [multipath antenna pathloss normalization]);
    
    realizationGroup = matlab.system.display.SectionGroup(...
        'Title', 'Realization', ...
        'Sections', randomization);
    
    groups = [mainGroup realizationGroup];
  end
end
end