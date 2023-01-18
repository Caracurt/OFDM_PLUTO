function [status] = tx_ofdm(IQmessage, center_frequency, sample_rate)
    %global plutoRadio
    tx_signal = IQmessage;

    %plutoRadio = findPlutoRadio;

    plutoTx = plutoTxInit(center_frequency,sample_rate);
    transmitRepeat(plutoTx,tx_signal);
    status = 1;
end