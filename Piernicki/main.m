%NADAJNIK
%Formatowanie danych do przesy³u
disp("Text to transmit: " + text);

textBinary = dec2bin(text, 8);
textLen = size(text,2);
bitsLen = 8*textLen;

bitsMatrix=zeros(textLen,8);
for m=1:textLen
    for k=1:8
        bitsMatrix(m,k)=str2double(textBinary(m,k));
    end
end

rowSum = sum(bitsMatrix,2);
rowSum = mod(rowSum', 2);
rowSum = rowSum.';
bitsMatrix(:,1) = rowSum;

disp("Text as binary matrix: ");
disp(bitsMatrix);

bitStream = bitsMatrix';
bitStream = bitStream(:)';

%Realizacja przesy³ania danych
% Amplitude for 0 bit
amplitude1 = 0; 
% Amplitude for 1 bit
amplitude2 = 1;

% Frequency of Modulating Signal
f = 2;
% Sampling rate - This will define the resoultion
fs = 100;
% Time for one bit
t = 0: 1/fs : bit_time;
% This time variable is just for plot
time = [];

ASK_signal = [];
Digital_signal = [];
for ii = 1: 1: length(bitStream)
    ASK_signal = [ASK_signal (bitStream(ii)==0)*amplitude1*sin(2*pi*f*t)+...
        (bitStream(ii)==1)*amplitude2*sin(2*pi*f*t)];
    
    time = [time t];
    t =  t + 1;
end

%Wykres sygna³u
axes(handles.axes1)
plot(ASK_signal,'LineWidth',1);
xlabel('n')
ylabel('zakodowane(n)')
title('Przebieg sygna³u w funkcji czasu');
grid  on;

%Obliczanie spektrum amplitudowego
spectrum = fft(ASK_signal);
spectrum_amp = abs(spectrum);

%Wykres spektrum aplitudowe
axes(handles.axes2)
plot(spectrum_amp);
xlabel('Spectral Index[Hz]');
ylabel('amplituda widma');
title('Widmo amplitudowe sygna³u');
grid on;

%DODANIE SZUMU
ASK_signal_with_noise = awgn(ASK_signal, snr);

%ODBIORNIK

%FILTRACJA
d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',0.15,'DesignMethod','butter');
ASK_signal_filtered = filtfilt(d1, ASK_signal_with_noise);

% Wykres sygna³u po filtracji szumu
axes(handles.axes3);
plot(time,ASK_signal_filtered,'LineWidth',1);
xlabel('Time[s]');
ylabel('Amplitude');
title('Przebieg sygna³u po filtracji w funkcji czasu');
grid  on;

bits_number = length(ASK_signal_filtered)/ length(t);
ASK_signal_splitted = reshape(ASK_signal_filtered, [length(t), bits_number]);
ASK_max_values = max(ASK_signal_splitted);
ASK_signal_decoded = round(ASK_max_values);

axes(handles.axes4);
scatter([1:length(ASK_signal_decoded)], ASK_signal_decoded);
xlabel('Bit');
ylabel('Wartoœæ');
title('Odebrany sygna³ bitowy');
grid  on;

ASK_signal_chars = reshape(ASK_signal_decoded, 8, length(ASK_signal_decoded)/8);
ASK_signal_chars = ASK_signal_chars.';

errorArray = mod(sum(ASK_signal_chars,2),2);
received_text = char(bin2dec(num2str(ASK_signal_chars))).';
errorsDetected = sum(errorArray);

% display values
handles.text8.String=received_text;
handles.text9.String=errorsDetected;