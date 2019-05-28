%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WLAN Link Modeling and performance analysis for different channel       %
% conditions                                                              %
%                                                                         %
% Author: Ata Niyazov (185112038)                                         %
%                                                                         %
% Work address: Kocaeli University                                        %
% Website: http://bilgisayar.kocaeli.edu.tr/                              %
% April 2019; Last revision: 04-Mart-2019                                 %
%                                                                         %
% Kocaeli University (C) Copyright 2019.                                  %
% All rights reserved.                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------- BEGIN CODE -------------------------------

%% END-TO-END simulation

%% Waveform Generation
% Create a format configuration object for a SISO VHT transmission
cfgVHT = wlanVHTConfig;            % Very High Throughput (VHT) designation
%VisualizeWLAN(cfgVHT);             % Visualize VHT Frame
cfgVHT.NumTransmitAntennas = 1;    % Transmit antennas
cfgVHT.NumSpaceTimeStreams = 1;    % Space-time streams
cfgVHT.APEPLength = 4096;          % APEP length in bytes
cfgVHT.MCS = 5;                    % Single spatial stream, 64-QAM
cfgVHT.ChannelBandwidth = 'CBW20'; % Transmitted signal bandwidth
Rs = wlanSampleRate(cfgVHT);       % Sampling rate

lstf = wlanLSTF(cfgVHT);  
lltf = wlanLLTF(cfgVHT);  
lsig = wlanLSIG(cfgVHT);

nonHTfield = [lstf;lltf;lsig]; % Combine the non-HT preamble fields

vhtsiga = wlanVHTSIGA(cfgVHT);
vhtstf = wlanVHTSTF(cfgVHT);
vhtltf = wlanVHTLTF(cfgVHT);
vhtsigb = wlanVHTSIGB(cfgVHT);

preamble = [lstf;lltf;lsig;vhtsiga;vhtstf;vhtltf;vhtsigb];

rng(0) % Initialize the random number generator
txPSDU = randi([0 1],cfgVHT.PSDULength*8,1); % Generate PSDU data in bits
data = wlanVHTData(txPSDU,cfgVHT);

txWaveform = [preamble;data]; % Transmit VHT PPDU

%% IEEE 802.11 Standard
% Standard: 802.11ac
% Bandwidth (MHz) 20, 40, 80, 160
% MIMO: Up to 8 spatial streams, MU-MIMO

%% Channel
% Parameterize the channel
tgacChannel = wlanTGacChannel;
tgacChannel.DelayProfile = 'Model-B';
tgacChannel.NumTransmitAntennas = cfgVHT.NumTransmitAntennas;
tgacChannel.NumReceiveAntennas = 1;
tgacChannel.LargeScaleFadingEffect = 'None';
tgacChannel.ChannelBandwidth = 'CBW20';
tgacChannel.TransmitReceiveDistance = 5;
tgacChannel.SampleRate = Rs;
tgacChannel.RandomStream = 'mt19937ar with seed';
tgacChannel.Seed = 10;

% Pass signal through the channel. Append zeroes to compensate for channel
% filter delay
txWaveform = [txWaveform;zeros(10,1)];
chanOut = tgacChannel(txWaveform);

snr = 40; % In dBs
rxWaveform = awgn(chanOut,snr,0);

%% Reciever
% Display the spectrum of the transmitted and received signals. The
% received signal spectrum is affected by the channel
spectrumAnalyzer  = dsp.SpectrumAnalyzer('SampleRate',Rs, ...
            'ShowLegend',true, ...
            'Window', 'Rectangular', ...
            'SpectralAverages',10, ...
            'YLimits',[-30 10], ... 
            'ChannelNames',{'Transmitted waveform','Received waveform'});
spectrumAnalyzer([txWaveform rxWaveform]);

%% Channel Estimation and Equalization
chInfo = info(tgacChannel); % Get characteristic information
% Channel filter delay, measured in samples 
chDelay  = chInfo.ChannelFilterDelay;
rxWaveform = rxWaveform(chDelay+1:end,:);

indField = wlanFieldIndices(cfgVHT);
indLLTF = indField.LLTF(1):indField.LLTF(2);
demodLLTF = wlanLLTFDemodulate(rxWaveform(indLLTF),cfgVHT);
% Estimate noise power in VHT fields
nVar = helperNoiseEstimate(demodLLTF,cfgVHT.ChannelBandwidth, ...
    cfgVHT.NumSpaceTimeStreams);

indVHTLTF = indField.VHTLTF(1):indField.VHTLTF(2);
demodVHTLTF = wlanVHTLTFDemodulate(rxWaveform(indVHTLTF,:),cfgVHT);
chanEstVHTLTF = wlanVHTLTFChannelEstimate(demodVHTLTF,cfgVHT);

figure
plot(20*log10(abs(chanEstVHTLTF)));
grid on;
title('Estimated Channel Response');
xlabel('Subcarrier index');
ylabel('Power (dB)');

indData = indField.VHTData(1):indField.VHTData(2);

% Recover the bits and equalized symbols in the VHT Data field using the
% channel estimates from VHT-LTF
[rxPSDU,~,eqSym] = wlanVHTDataRecover(rxWaveform(indData,:), ...
                    chanEstVHTLTF,nVar,cfgVHT);
        
% Compare transmit and receive PSDU bits       
numErr = biterr(txPSDU,rxPSDU);

% Plot equalized symbols
constellationDiagram = comm.ConstellationDiagram;
constellationDiagram.ReferenceConstellation = ...
    helperReferenceSymbols(cfgVHT);
% Compare received and reference constellation  
constellationDiagram(reshape(eqSym,[],1));      
constellationDiagram.Title = 'Equalized Data Symbols';

%------------------------------- END OF CODE ------------------------------

%% Analyze link performance by
% * Computing packet error rate
% * Bit error rate
% * Throughput measures
% ** Characterize and simulate wireless LAN fading path multichannel
% EVM & Spectral Emissions
% EVM BER PER