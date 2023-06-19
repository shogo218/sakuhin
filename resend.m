%%再送信するときの処理をする関数
function [rx,coarsePktOffset,pktOffset] = resend(txPSDU,cfg,AWGN,tgax,...
    mode,rx,coarsePktOffset,chanBW,ind,fs)

if mode==0
    tx = wlanWaveformGenerator(txPSDU,cfg);
    txPad = [tx; zeros(50,cfg.NumTransmitAntennas)];
    reset(tgax);
    rx = tgax(txPad);
    rx=AWGN(rx);
    coarsePktOffset = wlanPacketDetect(rx,chanBW);

elseif mode==1
    lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
    coarseFreqOff = wlanCoarseCFOEstimate(lstf,chanBW);
    rx = helperFrequencyOffset(rx,fs,-coarseFreqOff);

    nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
    finePktOffset = wlanSymbolTimingEstimate(nonhtfields,chanBW);

    pktOffset = coarsePktOffset+finePktOffset;
    
else
    error('Unexpected resend mode. Only 0 and 1 can use.')
end
end