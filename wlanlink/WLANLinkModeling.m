%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WLAN Link Modeling and performance analysis for different channel       %
% conditions                                                              %
%                                                                         %
% Author: Ata Niyazov (185112038), Fehime Yiğit (185112037)               %
%                                                                         %
% Work address: Kocaeli University                                        %
% Website: http://bilgisayar.kocaeli.edu.tr/                              %
% April 2019; Last revision: 21-May-2019                                  %
%                                                                         %
% Kocaeli University (C) Copyright 2019.                                  %
% All rights reserved.                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------- BEGIN CODE -------------------------------

close all; clear all; clc;

%% WLAN Link modellemesi ve farklı kanal koşulları için performans analizi

%% Giriş
% Projemizde MATLAB WLAN Toolbox kullanarak temel WLAN Link modelinin
% oluşturulması ve farklı koşullar altında performans analizini
% gerçekleştirilmesi göstermektedir. Bir IEEE 802.11ac VHT paketi
% oluşturduğumuz kanal üzerinden iletilir. Ardından alınan sinyal demodüle
% edilir ve iletilen bitleri geri kazanmak için eşitlenir (equalize) ve
% kodu çözülür (decode). Bu şekilde basit bir verici-kanal-alıcı
% simülasyonu gerçekleştirmektedir.

%% Dalga Biçimi Üretimi (Waveform Generation)
% Bu örnekte bir 802.11ac VHT aktarımı simüle edilmiştir. 802.11(TM)
% standardının VHT formatı için aktarım parametreleri, bir VHT
% konfigürasyon nesnesi kullanılarak konfigüre edilir ve "wlanVHTConfig"
% VHT yapılandırma nesnesi oluşturur. Bu örnekte, nesne 20 MHz kanal
% bant genişliği, MCS 5 ve tekli verici anten için yapılandırılmıştır.

% Create a format configuration object for a SISO VHT transmission
cfgVHT = wlanVHTConfig;            % Very High Throughput (VHT) designation
cfgVHT.NumTransmitAntennas = 1;    % Transmit antennas
%cfgVHT.NumTransmitAntennas = 8;    % 8 transmit antennas
cfgVHT.NumSpaceTimeStreams = 1;    % Space-time streams
%cfgVHT.NumSpaceTimeStreams = 8;    % 8 space-time streams
cfgVHT.APEPLength = 4096;          % APEP length in bytes
%cfgVHT.APEPLength = 3000;          % APEP length in bytes
% Modulation and coding scheme (MCS)
cfgVHT.MCS = 5;                    % Single spatial stream, 64-QAM
%cfgVHT.MCS = 9;                    % 256-QAM rate-5/6
% Channel bandwidth: 'CBW20','CBW40','CBW80','CBW160'
%cfgVHT.ChannelBandwidth = 'CBW20'; % Transmitted signal bandwidth
cfgVHT.ChannelBandwidth = 'CBW80'; % 80 MHz channel bandwidth
Rs = wlanSampleRate(cfgVHT);       % Sampling rate

%%
% Eğitim (training), sinyal (signal) ve veri (data) alanlarından oluşan
% tek bir VHT paketi üretilir:
%
% * Non-HT Short Training Field (L-STF)
% * Non-HT Long Training Field (L-LTF)
% * Non-HT Signal (L-SIG) field
% * VHT Signal A (VHT-SIG-A) field
% * VHT Short Training Field (VHT-STF)
% * VHT Long Training Field (VHT-LTF)
% * VHT Signal B (VHT-SIG-B) field
% * Data field
%
% Bu alanlar, WLAN Toolbox fonksiyonları kullanılarak ayrı ayrı
% üretilir ve bir VHT iletim paketi üretmek için birleştirilir.

%%
% Physical Protocol Data Unit (PPDU) ilk alanı L-STF’dir ve paket algılama
% (start of packet detection) ve otomatik kazanç kontrolü
% (AGC - automatic gain control) ayarı için kullanılır. Ayrıca ilk frekans
% ofset tahmini ve kaba zamanlama senkronizasyonu (time synchronization)
% için kullanılır. 'wlanLSTF' fonksiyonu, 'cfgVHT' yapılandırma nesnesinin
% içindeki parametrelerin bazılarını kullanarak zaman alanındaki
% (time-domain) L-STF alanını oluşturur.
lstf = wlanLSTF(cfgVHT);

%%
% L-LTF ince zaman senkronizasyonu (time synchronization), kanal kestirimi
% (channel estimation) ve ince frekans kayması kestirimi
% (frequency offset estimation) için kullanılır. 'wlanLLTF' fonksiyonu,
% zaman alanında (time-domain) L-LTF'yi oluşturur.
lltf = wlanLLTF(cfgVHT);

%%
% L-SIG alanı, non-HT formatı için veri hızı (data rate), modülasyon
% ve kod hızı (code rate) gibi paket konfigürasyonunu taşır. 'wlanLSIG'
% fonksiyonu, zaman alanındaki (time-domain) L-SIG alanını oluşturur.
lsig = wlanLSIG(cfgVHT);

%%
% L-STF, L-LTF ve L-SIG alanları VHT, HT-karışık (HT-mixed) ve HT-olmayan
% (non-HT) Orthogonal frequency-division multiplexing (OFDM) iletim
% biçimleri için ortaktır.
nonHTfield = [lstf;lltf;lsig]; % Combine the non-HT preamble fields

%%
% VHT için özgü sinyal ve eğitim alanları non-HT giriş alanlarından
% (preamble fields) sonra üretilir. VHT-SIG-A alanının amacı, alıcının
% veri yükü (data payload) kodunu çözmesine izin verecek bilgi sağlamaktır.
% VHT-SIG-A, VHT-SIG-A1 ve VHT-SIG-A2 sembollerinden oluşur. 'wlanVHTSIGA'
% fonksiyonu, zaman alanındaki (time-domain) VHT-SIG-A alanını oluşturur.
vhtsiga = wlanVHTSIGA(cfgVHT);

%%
% VHT-STF'nin amacı, bir MIMO aktarımındaki (multiple-input and
% multiple-output transmission) kazanç kontrolü tahminini (gain control
% estimation) iyileştirmek ve alıcının L-STF alanına benzer tekrar eden
% modeli tespit etmesine (detect the repeating pattern) yardımcı olmaktır.
% 'wlanVHTSTF' fonksiyonu, zaman alanındaki (time-domain) VHT-STF alanını
% oluşturur.
vhtstf = wlanVHTSTF(cfgVHT);

%%
% VHT-LTF, alıcı için verici ile alıcı arasındaki kanalı tahmin etmesi
% için bir araç sağlar. Uzay zaman akışlarının sayısına bağlı olarak
% (number of space time streams), 1,2,4,6 veya 8 VHT-LTF sembolünden
% oluşur. 'wlanVHTLTF' fonksiyonu, zaman alanındaki (time-domain)
% VHT-LTF'yi oluşturur.
vhtltf = wlanVHTLTF(cfgVHT);

%%
% VHT-SIG-B alanı, veri hızını (data rate) ve iletilen paketin veri alanı
% veri yükünün uzunluğunu (length of the data field payload) ayarlamak
% için kullanılır. 'wlanVHTSIGB' fonksiyonu, zaman alanındaki (time-domain)
% VHT-SIG-B alanını oluşturur.
vhtsigb = wlanVHTSIGB(cfgVHT);

%%
% Başlangıç (preamble) bölümünü, VHT formatı için oluşturulan sinyal ve
% eğitim alanlarıyla oluşturun.
preamble = [lstf;lltf;lsig;vhtsiga;vhtstf;vhtltf;vhtsigb];
plotVHTWaveform(preamble,cfgVHT,'Beamformed Preamble Wave with Fields Highlighted');

%%
% 'wlanVHTData' fonksiyonu, zaman alanında (time-domain) VHT veri alanını
% oluşturur. VHT formatı yapılandırması 'cfgVHT', veri alanını
% Physical Layer Convergence Procedure (PLCP) Service Data Unit (PSDU)
% bitlerinden üretmek için parametreleri belirtir. 'CfgVHT.PSDULength'
% özelliği, VHT veri alanında iletilecek bayt sayısını verir. Bu özellik,
% rasgele PSDU bit 'txPSDU' bit üretmek için kullanılır.

rng(0) % Initialize the random number generator
txPSDU = randi([0 1],cfgVHT.PSDULength*8,1); % Generate PSDU data in bits
data = wlanVHTData(txPSDU,cfgVHT);

% non-HT ve VHT başlangıç alanlarının (preamble) ardından veri alanını
% (data field) ekleyerek bir VHT dalga formu oluşturulur.
txWaveform = [preamble;data]; % Transmit VHT PPDU

%% VHT dalga biçimini alanları renklendirmiş olarak görüntüler
plotVHTWaveform(txWaveform,cfgVHT,'Beamformed TX Wave with Fields Highlighted');

%%
% Alternatif olarak, belirli bir format yapılandırması için dalga formu,
% tek bir fonksiyon çağrısı 'wlanWaveformGenerator' fonksiyonu kullanılarak
% da oluşturulabilir. Bu işlev bir veya daha fazla VHT paketi üretebilir.
% Varsayılan olarak OFDM pencereleme (OFDM windowing) oluşturulan dalga
% formuna (waveform) uygulanır. OFDM pencerelemesi hakkında daha fazla
% bilgi için 'wlanWaveformGenerator' belgelerine bakın.

%% Kanal bozuklukları (Channel Impairments)
% Bu bölüm kablosuz iletimin (over-the-air) etkilerini simüle eder.
% İletilen sinyal, kanal ve Additive white Gaussian noise (AWGN)
% tarafından bozulmuştur. AWGN'nin seviyesi dBs cinsinden verilmiştir.
% Bu örnekte, TGac kanal modeli, Model-B gecikme profili ile
% kullanılmıştır. Verici ve alıcı arasındaki mesafe 5 metreye eşit veya
% daha büyük olduğunda bu gecikme profili için, model Non-Line-of-Sight
% (N-LOS) konfigürasyonundadır. 'wlanTGacChannel' yardımında daha ayrıntılı
% olarak açıklanmaktadır.

% Parameterize the channel
tgacChannel = wlanTGacChannel;
%tgacChannel.DelayProfile = 'Model-B';
tgacChannel.DelayProfile = 'Model-D';
tgacChannel.NumTransmitAntennas = cfgVHT.NumTransmitAntennas;
tgacChannel.NumReceiveAntennas = 1;
%tgacChannel.NumReceiveAntennas = 8;
tgacChannel.LargeScaleFadingEffect = 'None';
%tgacChannel.ChannelBandwidth = 'CBW20';
tgacChannel.ChannelBandwidth = cfgVHT.ChannelBandwidth;
% Distance between the transmitter and receiver in meters, specified as
% a real positive scalar.
%tgacChannel.TransmitReceiveDistance = 5;
tgacChannel.TransmitReceiveDistance = 10; % Distance in meters for NLOS
tgacChannel.SampleRate = Rs;
tgacChannel.RandomStream = 'mt19937ar with seed';
tgacChannel.Seed = 10;
% Large-scale fading effects applied in the channel, specified as
% 'None', 'Pathloss', 'Shadowing', or 'Pathloss and shadowing'.
tgacChannel.LargeScaleFadingEffect = 'None';
%tgacChannel.LargeScaleFadingEffect = 'Pathloss';

% Sinyal kanaldan geçirilir. Kanal filtresi gecikmesini telafi etmek için
% elimizdeki dalga formunun (waveform) ardına sıfır eklenir.
txWaveform = [txWaveform;zeros(10,1)];
chanOut = tgacChannel(txWaveform);

snr = 30; % In dBs
rxWaveform = awgn(chanOut,snr,0);

%%
% Alınan sinyal spektrumu kanaldan etkilenir. İletilen ve alınan
% sinyallerin spektrumunu görüntülemek için:
spectrumAnalyzer  = dsp.SpectrumAnalyzer('SampleRate',Rs, ...
            'ShowLegend',true, ...
            'Window', 'Rectangular', ...
            'SpectralAverages',10, ...
            'YLimits',[-30 10], ... 
            'ChannelNames',{'Transmitted waveform','Received waveform'});
spectrumAnalyzer([txWaveform rxWaveform]);

%% Kanal Tahmini ve Denkleştirme (Channel Estimation and Equalization)
% Bu bölümde, zaman alanı (time-domain) VHT-LTF, alınan dalga formundan
% (received waveform) çıkarılır. Dalga formunun, kanal filtre gecikmesi
% hesaba katılarak paketin başlangıcına senkronize edildiği
% varsayılmaktadır. VHT-LTF demodüle edilir ve kanalı tahmin etmek için
% kullanılır. Alınan sinyal daha sonra VHT-LTF'den elde edilen kanal
% tahmini (channel estimate) kullanılarak eşitlenir (equalize).

%% TGac kanalı için 802.11ac Paket Hata Oranı(Packet Error Rate) simulasyonu
%% Simülasyon Parametreleri
% Vektördeki her SNR noktası için, bir paket sayısı üretilir, bir kanaldan
% geçirilir ve paket hata oranını belirlemek için demodüle edilir.

snrarr = 0:10:40;

%%
% Her bir SNR noktasında test edilen paketlerin sayısı iki parametre
% tarafından kontrol edilir:
%
% 'maxNumErrors' Her SNR noktasında simüle edilen maksimum paket hatası
% sayısıdır. Paket hatalarının sayısı bu sınıra ulaştığında, bu SNR
% noktasındaki simülasyon tamamlanmıştır.
%
% 'maxNumPackets' Her SNR noktasında simüle edilmiş maksimum paket
% sayısıdır ve paket hata limitine ulaşılmazsa simülasyon uzunluğunu sınırlar.
%
% Bu örnekte seçilen sayılar çok kısa bir simülasyona için yapılmıştır.
% Anlamlı sonuçlar için sayıları artırmanızı öneririz.

maxNumErrors = 10;   % The maximum number of packet errors at an SNR point
maxNumPackets = 100; % Maximum number of packets at an SNR point

%%
% Set the remaining variables for the simulation.

% Get the baseband sampling rate
%fs = wlanSampleRate(cfgVHT);

% Get the number of occupied subcarriers in VHT fields and FFT length
[vhtData,vhtPilots] = helperSubcarrierIndices(cfgVHT,'VHT');
Nst_vht = numel(vhtData)+numel(vhtPilots);
Nfft = helperFFTLength(cfgVHT);     % FFT length

% Set the sampling rate of the channel
%tgacChannel.SampleRate = fs;

% Indices for accessing each field within the time-domain packet
ind = wlanFieldIndices(cfgVHT);

%% SNR Noktalarını İşleme
% Her bir SNR noktası için birkaç paket test edilir ve paket hata oranı
% hesaplanır.
%
% Her paket için aşağıdaki işlem adımları gerçekleşir:
%
% # Tek bir paket dalga formu oluşturmak için bir PSDU oluşturulur
% ve kodlanır.
% # Dalga biçimi, TGac kanal modelinin farklı bir gerçekleştirmesinden
% geçer.
% # AWGN is added to the received waveform to create the desired average
% SNR per subcarrier after OFDM demodulation.
% # OFDM demodülasyonundan sonra alt taşıyıcı başına istenen ortalama
% SNR'yi oluşturmak için AWGN alınan dalga formuna eklenir.
% 'comm.AWGNChannel' doğru SNR sağlayacak şekilde yapılandırılmış.
% Yapılandırma, kanal içindeki normalizasyon için alıcı anten sayısı ve
% kullanılmayan alt taşıyıcılardaki OFDM demodülasyonu sırasında kaldırılan
% gürültü enerjisi ile ilgilidir.
% # Paket algılanır.
% # Kaba taşıyıcı frekans kayması tahmin edilir ve düzeltilir.
% # Hassas zamanlama senkronizasyonu yapıldı. L-STF, L-LTF ve L-SIG
% örnekleri, L-STF'nin başlangıcında veya sonunda paket saptamaya izin
% vermek için ince zamanlama için sağlanmıştır.
% # İnce taşıyıcı frekans kayması tahmin edildi ve düzeltildi.
% # VHT-LTF, senkronize edilmiş alınan dalga formundan çıkarılır.
% VHT-LTF, OFDM'nin demodüle edilmiş ve kanal kestirimi gerçekleştirilmiştir.
% # VHT Veri alanı, senkronize edilmiş alınan dalga formundan çıkarılır.
% PSDU, çıkartılan alan ve kanal tahmini kullanılarak kurtarılır.
% 
% SNR noktalarının işlenmesini paralelleştirmek için bir 'parfor' döngüsü
% kullanılabilir, bu nedenle her SNR noktası için bir AWGN kanalı
% oluşturulur ve 'comm.AWGNChannel' ile yapılandırılır. Paralel
% hesaplamanın artan hızda kullanılmasını sağlamak için 'for' ifadesine
% yorum yapın ve aşağıdaki 'parfor' ifadesini yorumdan çıkarın.

S = numel(snrarr);
packetErrorRate = zeros(S,1);
%parfor i = 1:S % Use 'parfor' to speed up the simulation
for i = 1:S     % Use 'for' to debug the simulation
    % Set random substream index per iteration to ensure that each
    % iteration uses a repeatable set of random numbers
    stream = RandStream('combRecursive','Seed',0);
    stream.Substream = i;
    RandStream.setGlobalStream(stream);
    
    % Create an instance of the AWGN channel per SNR point simulated
    awgnChannel = comm.AWGNChannel;
    awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
    % Normalization
    awgnChannel.SignalPower = 1/tgacChannel.NumReceiveAntennas; 
    % Account for energy in nulls
    awgnChannel.SNR = snrarr(i)-10*log10(Nfft/Nst_vht); 

    % Loop to simulate multiple packets
    numPacketErrors = 0;
    numPkt = 1; % Index of packet transmitted
    while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
        % Generate a packet waveform
        txPSDU = randi([0 1],cfgVHT.PSDULength*8,1); % PSDULength in bytes
        tx = wlanWaveformGenerator(txPSDU,cfgVHT);
        
        % Add trailing zeros to allow for channel delay
        tx = [tx; zeros(50,cfgVHT.NumTransmitAntennas)]; %#ok<AGROW>

        % Pass the waveform through the fading channel model
        reset(tgacChannel); % Reset channel for different realization
        rx = tgacChannel(tx);

        % Add noise
        rx = awgnChannel(rx);

        % Packet detect and determine coarse packet offset
        coarsePktOffset = wlanPacketDetect(rx,cfgVHT.ChannelBandwidth);
        if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
            numPacketErrors = numPacketErrors+1;
            numPkt = numPkt+1;
            continue; % Go to next loop iteration
        end
        
        % Extract L-STF and perform coarse frequency offset correction
        lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:); 
        coarseFreqOff = wlanCoarseCFOEstimate(lstf,cfgVHT.ChannelBandwidth);
        rx = helperFrequencyOffset(rx,Rs,-coarseFreqOff);
        
        % Extract the non-HT fields and determine fine packet offset
        nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:); 
        finePktOffset = wlanSymbolTimingEstimate(nonhtfields,...
                    cfgVHT.ChannelBandwidth);
        
        % Determine final packet offset
        pktOffset = coarsePktOffset+finePktOffset;
        
        % If packet detected outwith the range of expected delays from the
        % channel modeling; packet error
        if pktOffset>50
            numPacketErrors = numPacketErrors+1;
            numPkt = numPkt+1;
            continue; % Go to next loop iteration
        end

        % Extract L-LTF and perform fine frequency offset correction
        lltf = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:); 
        fineFreqOff = wlanFineCFOEstimate(lltf,cfgVHT.ChannelBandwidth);
        rx = helperFrequencyOffset(rx,Rs,-fineFreqOff);
        
        % Extract VHT-LTF samples from the waveform, demodulate and perform
        % channel estimation
        vhtltf = rx(pktOffset+(ind.VHTLTF(1):ind.VHTLTF(2)),:);
        vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,cfgVHT);
        
        % Get single stream channel estimate
        chanEstSSPilots = vhtSingleStreamChannelEstimate(vhtltfDemod,cfgVHT);
        
        % Channel estimate
        chanEst = wlanVHTLTFChannelEstimate(vhtltfDemod,cfgVHT);
        
        % Extract VHT Data samples from the waveform
        vhtdata = rx(pktOffset+(ind.VHTData(1):ind.VHTData(2)),:);
        
        % Estimate the noise power in VHT data field
        nVarVHT = vhtNoiseEstimate(vhtdata,chanEstSSPilots,cfgVHT);
        
        % Recover the transmitted PSDU in VHT Data
        rxPSDU = wlanVHTDataRecover(vhtdata,chanEst,nVarVHT,cfgVHT);
        
        % Determine if any bits are in error, i.e. a packet error
        packetError = any(biterr(txPSDU,rxPSDU));
        numPacketErrors = numPacketErrors+packetError;
        numPkt = numPkt+1;
    end

    % Calculate packet error rate (PER) at SNR point
    packetErrorRate(i) = numPacketErrors/(numPkt-1);
    disp(['SNR ' num2str(snrarr(i)) ' completed after ' ...
        num2str(numPkt-1) ' packets, PER: ' ... 
        num2str(packetErrorRate(i))]);
end

%% Plot Packet Error Rate vs SNR Results

figure
semilogy(snrarr,packetErrorRate,'-ob');
grid on;
xlabel('SNR (dB)');
ylabel('PER');
title('802.11ac 80MHz, MCS9, Direct Mapping, Channel Model D-NLOS');

%% 
% Bu örnekte, alınan sinyal bilinen bir kanal filtre gecikmesini telafi
% ederek paketin başlangıcına senkronize edilir.
chInfo = info(tgacChannel); % Get characteristic information
% Channel filter delay, measured in samples 
chDelay  = chInfo.ChannelFilterDelay;
rxWaveform = rxWaveform(chDelay+1:end,:);

%%
% Senkronizasyondan sonra, alıcının ilgili alanları alınan paketten
% çıkarması gerekir. 'wlanFieldIndices' fonksiyonu, bir paketteki ilk örneğe
% göre tüm alanların başlangıç ve bitiş zamanı etki alanı (time-domain)
% örnek dizinlerini döndürmek için kullanılır. Bu endeksler daha sonraki
% işlemler için gerekli alanları çıkarmak için kullanılır.
indField = wlanFieldIndices(cfgVHT);

%%
% Alınan OFDM sembollerinde Min Mean Square Error (MMSE) eşitlemesi yapmak
% için OFDM demodülasyonundan sonra gürültü gücünün bir tahmini gereklidir.
% Bu örnekte, VHT alanlarındaki gürültü gücü demodüle edilmiş L-LTF
% sembolleri kullanılarak tahmin edilir. L-LTF, alınan dalga formundan
% çıkarılır ve 'wlanLLTFDemodulate' fonksiyonu kullanılarak
% demodüle edilir.
indLLTF = indField.LLTF(1):indField.LLTF(2);
demodLLTF = wlanLLTFDemodulate(rxWaveform(indLLTF),cfgVHT);
% Estimate noise power in VHT fields
nVar = helperNoiseEstimate(demodLLTF,cfgVHT.ChannelBandwidth, ...
    cfgVHT.NumSpaceTimeStreams);

%%
% VHT-LTF'yi alınan sinyalden çıkarmak için başlangıç ve bitiş indeksleri
% bir indeks vektörü oluşturmak için kullanılır.
indVHTLTF = indField.VHTLTF(1):indField.VHTLTF(2);

%%
% VHT-LTF, antenden alınan tüm uzay-zaman akışları (space-time streams)
% arasındaki kanalı tahmin etmek için kullanılır. VHT-LTF, alınan dalga
% formundan çıkarılır ve 'wlanVHTLTFDemodulate' fonksiyonu kullanılarak
% demodüle edilir.
demodVHTLTF = wlanVHTLTFDemodulate(rxWaveform(indVHTLTF,:),cfgVHT);

%%
% Kanal tahmini, uygulanan uzaysal haritalamanın ve vericideki çoklu anten
% yapılandırması için döngüsel kaymaların etkisini içerir.
% 'wlanVHTLTFChannelEstimate' fonksiyonu, antenlerden alınan tüm
% uzay-zaman yayınları (space-time streams) arasındaki tahmini kanalı
% döndürür.
chanEstVHTLTF = wlanVHTLTFChannelEstimate(demodVHTLTF,cfgVHT);

%%
% İletim sinyali, aşağıdaki şekilde kanal frekans tepkisinde gösterildiği
% gibi derin bir solma (deep fade) ile karşılaşır. Kanal solmalarının
% etkisi, daha önce gösterilen spektrum çiziminde (spectrum plot) de
% görülebilir.
figure
plot(20*log10(abs(chanEstVHTLTF)));
grid on;
title('Estimated Channel Response');
xlabel('Subcarrier index');
ylabel('Power (dB)');

%%
% Veri alanını (data field) alınan sinyalden çıkarmak için, veri alanının
% başlangıç ve bitiş indeksleri bir indeks vektörü oluşturmak için
% kullanılır.
indData = indField.VHTData(1):indField.VHTData(2);

% VHT-LTF'den kanal tahminlerini (channel estimates) kullanarak bitleri
% ve VHT Veri alanındaki eşitlenmiş sembolleri (equalized symbols) kurtarın
[rxPSDU,~,eqSym] = wlanVHTDataRecover(rxWaveform(indData,:), ...
                    chanEstVHTLTF,nVar,cfgVHT);
        
% Verileri karşılaştırın ve PSDU bitlerini alarak, bit hata sayısı
% hesaplanır
disp('Number of Errors:');
numErr = biterr(txPSDU,rxPSDU)

%%
% Aşağıdaki grafik referans zayıflama (constellation) ile
% karşılaştırıldığında 'wlanVHTDataRecover' fonksiyonunun çıkışındaki
% eşitlenmiş sembollerin zayıflamasını göstermektedir. Kanal gürültüsünü
% arttırmak farklı zayıflama noktalarını yaymaya başlamalıdır.

% Plot equalized symbols
constellationDiagram = comm.ConstellationDiagram;
constellationDiagram.ReferenceConstellation = ...
    helperReferenceSymbols(cfgVHT);
% Compare received and reference constellation  
constellationDiagram(reshape(eqSym,[],1));      
constellationDiagram.Title = 'Equalized Data Symbols';

%% Analyze link performance by
% * Computing packet error rate (PER)
% * Bit error rate (BER)
% * Throughput measures

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperReferenceSymbols.m') helperReferenceSymbols.m>
% * <matlab:edit('helperNoiseEstimate.m') helperNoiseEstimate.m>

%% Referanslar
% # IEEE Std 802.11ac(TM)-2013 IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications -
% Amendment 4: Enhancements for Very High Throughput for Operation in Bands
% below 6 GHz.
% # Breit, G., H. Sampath, S. Vermani, et al. TGac Channel Model Addendum.
% Version 12. IEEE 802.11-09/0308r12, March 2010.

%------------------------------ END OF CODE -------------------------------