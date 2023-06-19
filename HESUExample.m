%% 802.11ax Packet Error Rate Simulation for Single-User Format
%
% This example shows how to measure the packet error rate of an IEEE(R)
% 802.11ax(TM) high efficiency (HE) single user format link.

% Copyright 2017-2019 The MathWorks, Inc.

%% Introduction
% In this example, an end-to-end simulation is used to determine the packet
% error rate for an 802.11ax [ <#10 1> ] single user format link for a
% selection of SNR points. At each SNR point, multiple packets are
% transmitted through a noisy TGax indoor channel, demodulated and the
% PSDUs recovered. The PSDUs are compared to those transmitted to determine
% the packet error rate. The processing for each packet is summarized in
% the following diagram.
%
% <<../HESUExampleDiagram.png>>
%
%% Waveform Configuration
% An HE single user (SU) packet is a full-band transmission to a single
% user. The transmit parameters for the HE SU format are configured using
% an <docid:wlan_ref#mw_12ce61d2-4b46-4827-9281-e65b2413e9ba
% wlanHESUConfig> object. The properties of the object contain the
% configuration. In this example, the object is configured for a 20 MHz
% channel bandwidth, 2 transmit antennas, 2 space-time streams, no space
% time block coding and 16-QAM rate-1/2 (MCS 3).

cfgHE = wlanHESUConfig;
cfgHE.ChannelBandwidth = 'CBW20';  % Channel bandwidth
cfgHE.NumSpaceTimeStreams = 2;     % Number of space-time streams
cfgHE.NumTransmitAntennas = 2;     % Number of transmit antennas
cfgHE.APEPLength = 1500;            % Payload length in bytes
cfgHE.ExtendedRange = false;       % Do not use extended range format
cfgHE.Upper106ToneRU = false;      % Do not use upper 106 tone RU
cfgHE.PreHESpatialMapping = false; % Spatial mapping of pre-HE fields
cfgHE.GuardInterval = 0.8;         % Guard interval duration
cfgHE.HELTFType = 4;               % HE-LTF compression mode
cfgHE.ChannelCoding = 'BCC';      % Channel coding
cfgHE.MCS = 3;                     % Modulation and coding scheme

%提案手法を比較するためのもう一つの無線LANの設定
cfgHE2 = wlanHESUConfig;
cfgHE2.ChannelBandwidth = 'CBW20';  % Channel bandwidth
cfgHE2.NumSpaceTimeStreams = 2;     % Number of space-time streams
cfgHE2.NumTransmitAntennas = 2;     % Number of transmit antennas
cfgHE2.APEPLength = 1500;            % Payload length in bytes
cfgHE2.ExtendedRange = false;       % Do not use extended range format
cfgHE2.Upper106ToneRU = false;      % Do not use upper 106 tone RU
cfgHE2.PreHESpatialMapping = false; % Spatial mapping of pre-HE fields
cfgHE2.GuardInterval = 0.8;         % Guard interval duration
cfgHE2.HELTFType = 4;               % HE-LTF compression mode
cfgHE2.ChannelCoding = 'BCC';      % Channel coding
cfgHE2.MCS = 3;                     % Modulation and coding scheme

%% Channel Configuration
% In this example, a TGax NLOS indoor channel model is used with delay
% profile Model-B. Model-B is considered NLOS when the distance between
% transmitter and receiver is greater than or equal to 5 meters. This is
% described further in
% <docid:wlan_ref#mw_43b5900e-69e1-4636-b084-1e72dbd46293 wlanTGaxChannel>.
% A 2x2 MIMO channel is simulated in this example.

% Create and configure the TGax channel
chanBW = cfgHE.ChannelBandwidth;
tgaxChannel = wlanTGaxChannel;
tgaxChannel.DelayProfile = 'Model-B';
tgaxChannel.NumTransmitAntennas = cfgHE.NumTransmitAntennas;
tgaxChannel.NumReceiveAntennas = 2;
tgaxChannel.TransmitReceiveDistance = 5; % Distance in meters for NLOS
tgaxChannel.ChannelBandwidth = chanBW;
tgaxChannel.LargeScaleFadingEffect = 'None';
fs = wlanSampleRate(cfgHE);
tgaxChannel.SampleRate = fs;

%% Simulation Parameters
% For each SNR point (dB) in the |snr| vector a number of packets are
% generated, passed through a channel and demodulated to determine the
% packet error rate.

snr = 0:5:50;

%%
% The number of packets tested at each SNR point is controlled by two
% parameters:
%
% # |maxNumErrors| is the maximum number of packet errors simulated at each
% SNR point. When the number of packet errors reaches this limit, the
% simulation at this SNR point is complete.
% # |maxNumPackets| is the maximum number of packets simulated at each SNR
% point and limits the length of the simulation if the packet error limit
% is not reached.
%
% The numbers chosen in this example will lead to a very short simulation.
% For statistically meaningful results we recommend increasing these
% numbers.

%エラーの数関係なく全てのパケットを送信
maxNumPackets = 10; % The Maximum number of packets at an SNR point

%% Processing SNR Points
% For each SNR point a number of packets are tested and the packet error
% rate calculated. The pre-HE preamble of 802.11ax is backwards compatible
% with 802.11ac(TM), therefore in this example the front-end
% synchronization components for a VHT waveform are used to synchronize the
% HE waveform at the receiver. For each packet the following processing
% steps occur:
%
% # A PSDU is created and encoded to create a single packet waveform.
% # The waveform is passed through an indoor TGax channel model. Different
% channel realizations are modeled for different packets.
% # AWGN is added to the received waveform to create the desired average
% SNR per subcarrier after OFDM demodulation. The <docid:comm_ref#buiamu7-1
% comm.AWGNChannel> object is configured to provide the correct SNR. The
% configuration accounts for normalization within the channel by the number
% of receive antennas, and the noise energy in unused subcarriers which are
% removed during OFDM demodulation.
% # The packet is detected.
% # Coarse carrier frequency offset is estimated and corrected.
% # Fine timing synchronization is established. The L-STF, L-LTF and L-SIG
% samples are provided for fine timing to allow for packet detection at the
% start or end of the L-STF.
% # Fine carrier frequency offset is estimated and corrected.
% # The HE-LTF is extracted from the synchronized received waveform. The
% HE-LTF is OFDM demodulated and channel estimation is performed.
% # The data field is extracted from the synchronized received waveform
% and OFDM demodulated.
% # Common phase error pilot tracking is performed to track any residual
% carrier frequency offset.
% # Noise estimation is performed using the demodulated data field pilots
% and single-stream channel estimate at pilot subcarriers.
% # The phase corrected OFDM symbols are equalized with the channel
% estimate.
% # The equalized symbols are demodulated and decoded to recover the PSDU.
%
% A <docid:matlab_ref#f71-813245 parfor> loop can be used to parallelize
% processing of the SNR points, therefore for each SNR point an AWGN
% channel is created and configured by using the |comm.AWGNChannel| object.
% To enable the use of parallel computing for increased speed comment out
% the 'for' statement and uncomment the 'parfor' statement below.

numSNR = numel(snr); % Number of SNR points
packetErrorRate = zeros(1,numSNR);
throughput = zeros(1,numSNR);
fs = wlanSampleRate(cfgHE);

packetErrorRate_harq = zeros(1,numSNR);
throughput_harq = zeros(1,numSNR);

% Get occupied subcarrier indices and OFDM parameters
ofdmInfo = wlanHEOFDMInfo('HE-Data',cfgHE);

% Indices to extract fields from the PPDU
ind = wlanFieldIndices(cfgHE);

%parfor isnr = 1:numSNR % Use 'parfor' to speed up the simulation
for isnr = 1:numSNR
    % Set random substream index per iteration to ensure that each
    % iteration uses a repeatable set of random numbers
    %stream = RandStream('combRecursive','Seed',99);
    %stream.Substream = isnr;
    %RandStream.setGlobalStream(stream);

    % Create an instance of the AWGN channel per SNR point simulated
    awgnChannel = comm.AWGNChannel;
    awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
    awgnChannel.SignalPower = 1/tgaxChannel.NumReceiveAntennas;
    % Account for noise energy in nulls so the SNR is defined per
    % active subcarrier
    awgnChannel.SNR = snr(isnr)-10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);
    SNR = awgnChannel.SNR;

    % Loop to simulate multiple packets
    numPacketErrors = 0;
    numPkt = 1; % Index of packet transmitted
    numPacketErrors2 = 0;
    numPkt2 = 1;
    
    %IR-HARQの設定
    cfgHARQ=struct;
    cfgHARQ.Symbol = 'BPSK';
    cfgHARQ.rate = '5/6';
    cfgHARQ.trellis = 0;
    cfgHARQ.puncpat = 0;
    cfgHARQ.NACK = 0;
    cfgHARQ.Ndata = 0;
    
    %SNRの値によってディジタル変調を変更
    if SNR <= 14 && 8 < SNR
        cfgHARQ.Symbol = 'QPSK';
    elseif 14 < SNR && SNR <= 20
        cfgHARQ.Symbol = '16QAM';
    elseif 20<SNR
        cfgHARQ.Symbol = '64QAM';
    end
    

    S_sum  = 0;
    maxRetrans = 3; %最大再送回数
    while numPkt <= maxNumPackets
        %送信失敗した場合の状態を管理
        errorFlag = 0;
        errorFlag2 = 0;

        % Generate a packet with random PSDU
        psduLength = getPSDULength(cfgHE); % PSDU length in bytes
        txPSDU = randi([0 1],psduLength*8,1);

        [txPSDU_harq,trellis,puncpat] = harq_launcher(txPSDU,cfgHARQ,...
            1,psduLength*8);
        cfgHARQ.trellis = trellis;
        cfgHARQ.puncpat = puncpat;

        %IR-HARQによってデータサイズが変更されたためバッファデータ追加による調整
        blankSize = mod(length(txPSDU_harq),8);
        if blankSize ~= 0
            blankSize = 8 - blankSize;
        end
        harqLength = length(txPSDU_harq) + blankSize;
        cfgHE2.APEPLength = (harqLength)/8;
        psduLength2 = getPSDULength(cfgHE2);
        blank = zeros((psduLength2*8-length(txPSDU_harq)),1);
        txPSDU_harq = cat(1,txPSDU_harq,blank);
        ind2 = wlanFieldIndices(cfgHE2);
        
        tx = wlanWaveformGenerator(txPSDU,cfgHE);
        tx_harq = wlanWaveformGenerator(txPSDU_harq,cfgHE2);
        
        % Add trailing zeros to allow for channel delay
        txPad = [tx; zeros(50,cfgHE.NumTransmitAntennas)];
        txPad_harq = [tx_harq; zeros(50,cfgHE2.NumTransmitAntennas)];

        % Pass through a fading indoor TGax channel
        reset(tgaxChannel); % Reset channel for different realization
        rx = tgaxChannel(txPad);
        reset(tgaxChannel);
        rx_harq = tgaxChannel(txPad_harq);

        % Pass the waveform through AWGN channel
        rx = awgnChannel(rx);
        rx_harq = awgnChannel(rx_harq);

        % Packet detect and determine coarse packet offset
        coarsePktOffset = wlanPacketDetect(rx,chanBW);
        coarsePktOffset2 = wlanPacketDetect(rx_harq,chanBW);
        if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
            numPacketErrors = numPacketErrors + 1;
            numPkt = numPkt + 1;
            errorFlag = 1; % Go to next loop iteration
        end
        if isempty(coarsePktOffset2)
            %再送信処理
            errorFlag2 = errorFlag2+1;
            if errorFlag2 <= maxRetrans
                for i = errorFlag2:maxRetrans
                    %再送信した回数によって符号化率を変更
                    if errorFlag2 == 1
                        cfgHARQ.rate = 3/4;
                    elseif errorFlag2 == 2
                        cfgHARQ.rate = 2/3;
                    else
                        cfgHARQ.rate = 1/2;
                    end

                    [txPSDU_harq,~,puncpat] = harq_launcher(txPSDU,cfgHARQ,...
                        1,psduLength*8);

                    %IR-HARQによってデータサイズが変更されたためバッファデータ追加による調整
                    blankSize = mod(length(txPSDU_harq),8);
                    if blankSize ~= 0
                        blankSize = 8 - blankSize;
                    end
                    harqLength = length(txPSDU_harq) + blankSize;
                    cfgHE2.APEPLength = (harqLength) / 8;
                    psduLength2 = getPSDULength(cfgHE2);
                    blank = zeros((psduLength2*8-length(txPSDU_harq)),1);
                    txPSDU_harq = cat(1,txPSDU_harq,blank);

                    ind2 = wlanFieldIndices(cfgHE2);
                    [rx_harq,coarsePktOffset2] = resend(txPSDU_harq,cfgHE2,awgnChannel,...
                        tgaxChannel,0,0,0,chanBW,ind2,fs);

                    if isempty(coarsePktOffset2)
                       errorFlag2 = errorFlag2 + 1;
                       if errorFlag2 >= maxRetrans
                           packetError2 = 1;
                           numPacketErrors2 = numPacketErrors2 + 1;
                           numPkt2 = numPkt2 + 1;
                       end
                    else
                        break;
                    end
                end
                cfgHARQ.puncpat = puncpat;
            else
                numPacketErrors2 = numPacketErrors2 + 1;
                numPkt2 = numPkt2 + 1;
            end
        end
        
        %従来手法の処理
        if errorFlag == 0
            % Extract L-STF and perform coarse frequency offset correction
            lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
            coarseFreqOff = wlanCoarseCFOEstimate(lstf,chanBW);
            rx = helperFrequencyOffset(rx,fs,-coarseFreqOff);

            % Extract the non-HT fields and determine fine packet offset
            nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
            finePktOffset = wlanSymbolTimingEstimate(nonhtfields,chanBW);

            % Determine final packet offset
            pktOffset = coarsePktOffset + finePktOffset;
            % If packet detected outwith the range of expected delays from
            % the channel modeling; packet error
            if pktOffset > 50
                numPacketErrors = numPacketErrors + 1;
                numPkt = numPkt + 1;
                errorFlag = 1; % Go to next loop iteration
            end
        end
        
        %提案手法の処理
        if errorFlag2 <= maxRetrans
            % Extract L-STF and perform coarse frequency offset correction
            lstf_harq = rx_harq(coarsePktOffset2+(ind2.LSTF(1):ind2.LSTF(2)),:);
            coarseFreqOff2 = wlanCoarseCFOEstimate(lstf_harq,chanBW);
            rx_harq = helperFrequencyOffset(rx_harq,fs,-coarseFreqOff2);

            % Extract the non-HT fields and determine fine packet offset
            nonhtfields_harq = rx_harq(coarsePktOffset2+(ind2.LSTF(1):ind2.LSIG(2)),:);
            finePktOffset2 = wlanSymbolTimingEstimate(nonhtfields_harq,chanBW);

            % Determine final packet offset
            pktOffset_harq = coarsePktOffset2 + finePktOffset2;
            if pktOffset_harq > 50
                %再送信処理
                errorFlag2 = errorFlag2+1;
                if errorFlag2 <= maxRetrans
                    for i = errorFlag2:maxRetrans
                        %再送信した回数によって符号化率を変更
                        if errorFlag2 == 1
                            cfgHARQ.rate = 3/4;
                        elseif errorFlag2 == 2
                            cfgHARQ.rate = 2/3;
                        else
                            cfgHARQ.rate = 1/2;
                        end
                        [txPSDU_harq,~,puncpat] = harq_launcher(txPSDU,cfgHARQ,...
                            1,psduLength*8);

                        %IR-HARQによってデータサイズが変更されたためバッファデータ追加による調整
                        blankSize = mod(length(txPSDU_harq),8);
                        if blankSize ~= 0
                            blankSize = 8 - blankSize;
                        end
                        harqLength = length(txPSDU_harq) + blankSize;
                        cfgHE2.APEPLength = (harqLength)/8;
                        psduLength2 = getPSDULength(cfgHE2);
                        blank = zeros((psduLength2*8-length(txPSDU_harq)),1);
                        txPSDU_harq = cat(1,txPSDU_harq,blank);
                        
                        ind2 = wlanFieldIndices(cfgHE2);
                        [rx_harq,coarsePktOffset2] = resend(txPSDU_harq,cfgHE2,awgnChannel,...
                            tgaxChannel,0,0,0,chanBW,ind2,fs);
                        if isempty(coarsePktOffset2)
                            errorFlag2 = errorFlag2 + 1;
                            continue;
                        else
                            [rx_harq,~,pktOffset_harq] = resend(rx_harq,cfgHE2,awgnChannel,...
                                tgaxChannel,1,rx_harq,coarsePktOffset2,chanBW,ind2,fs);
                            if pktOffset_harq > 50
                                errorFlag2 = errorFlag2 + 1;
                                if errorFlag2 >= maxRetrans
                                    packetError2 = 1;
                                    numPacketErrors2 = numPacketErrors2 + 1;
                                    numPkt2 = numPkt2+1;
                                end
                                continue;
                            else
                                break;
                            end
                        end
                    end
                    cfgHARQ.puncpat = puncpat;
                else
                    packetError2 = 1;
                    numPacketErrors2 = numPacketErrors2 + 1;
                    numPkt2 = numPkt2 + 1;
                end
            end
        end

        %従来手法の処理
        if errorFlag == 0
            % Extract L-LTF and perform fine frequency offset correction
            rxLLTF = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
            fineFreqOff = wlanFineCFOEstimate(rxLLTF,chanBW);
            rx = helperFrequencyOffset(rx,fs,-fineFreqOff);

            % HE-LTF demodulation and channel estimation
            rxHELTF = rx(pktOffset+(ind.HELTF(1):ind.HELTF(2)),:);
            heltfDemod = wlanHEDemodulate(rxHELTF,'HE-LTF',cfgHE);
            [chanEst,pilotEst] = heLTFChannelEstimate(heltfDemod,cfgHE);

            % Data demodulate
            rxData = rx(pktOffset+(ind.HEData(1):ind.HEData(2)),:);
            demodSym = wlanHEDemodulate(rxData,'HE-Data',cfgHE);

            % Pilot phase tracking
            demodSym = heCommonPhaseErrorTracking(demodSym,chanEst,cfgHE);

            % Estimate noise power in HE fields
            nVarEst = heNoiseEstimate(demodSym(ofdmInfo.PilotIndices,:,:),pilotEst,cfgHE);

            % Extract data subcarriers from demodulated symbols and channel
            % estimate
            demodDataSym = demodSym(ofdmInfo.DataIndices,:,:);
            chanEstData = chanEst(ofdmInfo.DataIndices,:,:);

            % Equalization and STBC combining
            [eqDataSym,csi] = heEqualizeCombine(demodDataSym,chanEstData,nVarEst,cfgHE);

            % Recover data
            rxPSDU = wlanHEDataBitRecover(eqDataSym,nVarEst,csi,cfgHE);

            % Determine if any bits are in error, i.e. a packet error
            packetError = ~isequal(txPSDU,rxPSDU);
            numPacketErrors = numPacketErrors+packetError;
            numPkt = numPkt+1;
        end
        
        %提案手法の処理
        if errorFlag2 <= maxRetrans
            % Extract L-LTF and perform fine frequency offset correction
            rxLLTF_harq = rx_harq(pktOffset_harq+(ind2.LLTF(1):ind2.LLTF(2)),:);
            fineFreqOff2 = wlanFineCFOEstimate(rxLLTF_harq,chanBW);
            rx_harq = helperFrequencyOffset(rx_harq,fs,-fineFreqOff2);

            % HE-LTF demodulation and channel estimation
            rxHELTF_harq = rx_harq(pktOffset_harq+(ind2.HELTF(1):ind2.HELTF(2)),:);
            heltfDemod_harq = wlanHEDemodulate(rxHELTF_harq,'HE-LTF',cfgHE2);
            [chanEst2,pilotEst2] = heLTFChannelEstimate(heltfDemod_harq,cfgHE2);

            % Data demodulate
            rxData_harq = rx_harq(pktOffset_harq+(ind2.HEData(1):ind2.HEData(2)),:);
            demodSym_harq = wlanHEDemodulate(rxData_harq,'HE-Data',cfgHE2);

            % Pilot phase tracking
            demodSym_harq = heCommonPhaseErrorTracking(demodSym_harq,chanEst2,cfgHE2);

            % Estimate noise power in HE fields
            nVarEst_harq = heNoiseEstimate(demodSym_harq(ofdmInfo.PilotIndices,:,:),pilotEst2,cfgHE2);

            % Extract data subcarriers from demodulated symbols and channel
            % estimate
            demodDataSym_harq = demodSym_harq(ofdmInfo.DataIndices,:,:);
            chanEstData_harq = chanEst2(ofdmInfo.DataIndices,:,:);

            % Equalization and STBC combining
            [eqDataSym_harq,csi2] = heEqualizeCombine(demodDataSym_harq,...
                chanEstData_harq,nVarEst_harq,cfgHE2);

            % Recover data
            rxPSDU_harq = wlanHEDataBitRecover(eqDataSym_harq,nVarEst_harq,...
                csi2,cfgHE2);
            rxPSDU_harq1 = rxPSDU_harq(1:(harqLength-blankSize),1);
            rxPSDU_harq=harq_launcher(double(rxPSDU_harq1),cfgHARQ,...
            2,psduLength*8);
            
            % Determine if any bits are in error, i.e. a packet error
            packetError2 = ~isequal(txPSDU,rxPSDU_harq);

            %復号失敗した場合の処理
            if packetError2 == 1
                cfgHARQ.NACK = 1;
                if errorFlag2 <= maxRetrans
                    for i = errorFlag2:maxRetrans

                        %再送信した回数によって符号化率を変更
                        if errorFlag2 == 1
                            cfgHARQ.rate = 3/4;
                        elseif errorFlag2 == 2
                            cfgHARQ.rate = 2/3;
                        else
                            cfgHARQ.rate = 1/2;
                        end
                        [txPSDU_harq,~,puncpat] = harq_launcher(txPSDU,...
                            cfgHARQ,1,psduLength*8);

                        %IR-HARQによってデータサイズが変更されたためバッファデータ追加による調整
                        blankSize = mod(length(txPSDU_harq),8);
                        if blankSize ~= 0
                            blankSize = 8 - blankSize;
                        end
                        harqLength = length(txPSDU_harq) + blankSize;
                        cfgHE2.APEPLength = (harqLength) / 8;
                        psduLength2 = getPSDULength(cfgHE2);
                        blank = zeros((psduLength2*8-length(txPSDU_harq)),1);
                        txPSDU_harq = cat(1,txPSDU_harq,blank);
                        
                        ind2 = wlanFieldIndices(cfgHE2);
                        cfgHARQ.puncpat = puncpat;
                        [rx_harq,coarsePktOffset2] = resend(txPSDU_harq, cfgHE2,...
                            awgnChannel, tgaxChannel, 0, 0, 0, chanBW, ind2, fs);
                        if isempty(coarsePktOffset2)
                            errorFlag2 = errorFlag2 + 1;
                            continue;
                        end
                        [rx_harq,~,pktOffset_harq] = resend(rx_harq, cfgHE2, awgnChannel,...
                            tgaxChannel, 1, rx_harq, coarsePktOffset2, chanBW, ind2, fs);
                        if pktOffset_harq > 50
                            errorFlag2 = errorFlag2 + 1;
                            continue;
                        end
                        % Extract L-LTF and perform fine frequency offset correction
                        rxLLTF_harq = rx_harq(pktOffset_harq+(ind2.LLTF(1):ind2.LLTF(2)),:);
                        fineFreqOff2 = wlanFineCFOEstimate(rxLLTF_harq,chanBW);
                        rx_harq = helperFrequencyOffset(rx_harq,fs,-fineFreqOff2);

                        % HE-LTF demodulation and channel estimation
                        rxHELTF_harq = rx_harq(pktOffset_harq+(ind2.HELTF(1):ind2.HELTF(2)),:);
                        heltfDemod_harq = wlanHEDemodulate(rxHELTF_harq,'HE-LTF',cfgHE2);
                        [chanEst2,pilotEst2] = heLTFChannelEstimate(heltfDemod_harq,cfgHE2);

                        % Data demodulate
                        rxData_harq = rx_harq(pktOffset_harq+(ind2.HEData(1):ind2.HEData(2)),:);
                        demodSym_harq = wlanHEDemodulate(rxData_harq,'HE-Data',cfgHE2);

                        % Pilot phase tracking
                        demodSym_harq = heCommonPhaseErrorTracking(demodSym_harq,chanEst2,cfgHE2);

                        % Estimate noise power in HE fields
                        nVarEst_harq = heNoiseEstimate(demodSym_harq(ofdmInfo.PilotIndices,:,:),pilotEst2,cfgHE2);

                        % Extract data subcarriers from demodulated symbols and channel
                        % estimate
                        demodDataSym_harq = demodSym_harq(ofdmInfo.DataIndices,:,:);
                        chanEstData_harq = chanEst2(ofdmInfo.DataIndices,:,:);

                        % Equalization and STBC combining
                        [eqDataSym_harq,csi2] = heEqualizeCombine(demodDataSym_harq,...
                            chanEstData_harq,nVarEst_harq,cfgHE2);
                        
                        % Recover data
                        rxPSDU_harq = wlanHEDataBitRecover(eqDataSym_harq,nVarEst_harq,...
                            csi2,cfgHE2);
                        rxPSDU_harq1 = rxPSDU_harq(1:(harqLength-blankSize),1);
                        rxPSDU_harq = harq_launcher(rxPSDU_harq1,cfgHARQ,...
                            2,psduLength*8);
                        packetError2 = ~isequal(txPSDU,rxPSDU_harq);
                        errorFlag2 = errorFlag2+packetError2;
                        if packetError2 == 0
                            numPacketErrors2 = numPacketErrors2 + packetError2;
                            numPkt2 = numPkt2 + 1;
                            break;
                        end
                    end
                    if packetError2 == 1
                        numPacketErrors2 = numPacketErrors2 + packetError2;
                        numPkt2 = numPkt2 + 1;
                    end
                else
                    numPacketErrors2 = numPacketErrors2 + packetError2;
                    numPkt2 = numPkt2 + 1;
                end
            else
                numPacketErrors2 = numPacketErrors2 + packetError2;
                numPkt2 = numPkt2 + 1;
            end
        end
        txtime2 = (length(tx_harq)/fs) * 1e6;
        S = ((1 - packetError2) * length(txPSDU)) / (txtime2);
        S_sum = S_sum + S;
    end
    succ = (numPkt - 1) - numPacketErrors;
    succ2 = (numPkt2 - 1) - numPacketErrors2;
    txtime = (length(tx)/fs) * 1e6;
    %スループット(通信速度)の計算
    throughput(isnr) = 1e-6 * ((maxNumPackets - numPacketErrors)...
        *(length(txPSDU))/(maxNumPackets * (txtime) * 1e-6));
    throughput_harq(isnr) = 2 * S_sum / (numPkt2-1);
    % Calculate packet error rate (PER) at SNR point
    packetErrorRate(isnr) = numPacketErrors / (numPkt - 1);
    packetErrorRate_harq(isnr) = numPacketErrors2 / (numPkt2 - 1);
    disp(['MCS ' num2str(cfgHE.MCS) ','...
          ' SNR ' num2str(snr(isnr)) ...
          ' completed after ' num2str(numPkt - 1) ' packets,'...
          ' PER:' num2str(packetErrorRate(isnr)) ' ,Throuput '...
          num2str(throughput(isnr)),' [Mbps] ',num2str(succ)]);
    disp('Used IR-HARQ');
    disp(['MCS ' num2str(cfgHE.MCS) ','...
          ' SNR ' num2str(snr(isnr)) ...
          ' completed after ' num2str(numPkt - 1) ' packets,'...
          ' PER:' num2str(packetErrorRate_harq(isnr)) ' ,Throuput '...
          num2str(throughput_harq(isnr)),' [Mbps] ',num2str(succ2)]);
end

%% Plot Packet Error Rate vs SNR

figure;
semilogy(snr,packetErrorRate,'-b^',snr,packetErrorRate_harq,'-ro','LineWidth',1);
hold on;
grid on;
xlabel('SNR (dB)');
ylabel('PER (Packet Error Rate)');
xlim([0 50])
ylim([1.e-01 1.e+00])
legend('Conventional','Proposed');

figure;
plot(snr,throughput,'-b^',snr,throughput_harq,'-ro','LineWidth',1);
hold on;
grid on;
xlabel('SNR (dB)');
ylabel('Throughput (Mbps)');
xlim([0 50])
ylim([0 150])
legend('Conventional','Proposed');

%%
% The number of packets tested at each SNR point is controlled by two
% parameters: |maxNumErrors| and |maxNumPackets|. For meaningful results,
% these values should be larger than those presented in this example. As an
% example, the figure below was created by running a longer simulation with
% |maxNumErrors|:1e3 and |maxNumPackets|:1e4.
%
% <<../HESUExamplePER.png>>

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('heCommonPhaseErrorTracking.m') heCommonPhaseErrorTracking.m>
% * <matlab:edit('heEqualizeCombine.m') heEqualizeCombine.m>
% * <matlab:edit('helperFrequencyOffset.m') helperFrequencyOffset.m>
% * <matlab:edit('heLTFChannelEstimate.m') heLTFChannelEstimate.m>
% * <matlab:edit('heNoiseEstimate.m') heNoiseEstimate.m>

%% Selected Bibliography
% # IEEE P802.11ax(TM)/D4.1 Draft Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications -
% Amendment 6: Enhancements for High Efficiency WLAN.
