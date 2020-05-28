function output = preprocessSuitWidar(csi,s)
%% Input validation
%   - The parameters have been validated in preprocess()
%   - Only validating the suit name whether it is 'Widar'
%
%   Input description:
%   csi - a complex matrix with size of [packetNum, subcarrierNum*links]
%   s - a struct of processing suite
assert(contains(s.suite,'Widar','IgnoreCase',true),'The suite of preprocessSuitWidar() should be Widar');

%% Noise Filtering using Butterworth Filter
if contains(s.filter,'bpf')
    csiFilter = BPF(csiCM,s.fc,s.fs); % Band-pass filter
elseif contains(s.filter,'lpf')
    csiFilter = LPF(csiCM,s.fc,s.fs); % Low-pass filter
elseif contains(s.filter,'hpf')
    csiFilter = HPF(csiCM,s.fc,s.fs); % High-pass filter
end

%% PCA-based denoising
csiPCA = PCA(csiFilter,s.pca);

%% Spectrogram Generation
[output.psd,output.frequency,output.time] = Spectrogram_Generation( ...
    csiPCA, s.stft(1), s.stft(1)-s.stft(2), s.fs, s.stft(3));

%% PSD coutour extraction
[row, column] = size(output.psd);
    
threshhold1 = quantile(output.psd(:)',0.75);
contour = zeros(1,column);
for i=1:column
    couttourPoint = find(output.psd(:,i)>=threshhold1,1,'last');
    if isempty(couttourPoint)
        coutour(1,i) = 0;
    else
        coutour(1,i) = couttourPoint;
    end
end
%smooth the coutour
halfWin=floor(s.winSize/2);
coutour = [zeros(1,halfWin), coutour, zeros(1,halfWin)];
for i= halfWin + 1:column+halfWin
    smoothCoutour(1,i-halfWin) = sum(coutour(1,i-halfWin:i+halfWin))/(halfWin*2+1);
end
output.coutour = smoothCoutour;
end

