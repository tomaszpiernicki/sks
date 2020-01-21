%text = 'Test';

%NADAJNIK
%Formatowanie danych do przesy�u
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

%Realizacja przesy�ania danych
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

%Wykres sygna�u
% subplot(5,1,1);
axes(handles.axes1)
plot(ASK_signal,'LineWidth',1);
xlabel('n')
ylabel('zakodowane(n)')
title('Przebieg sygna�u w funkcji czasu');
%axis([0 time(end) 1.5 1.5]);
grid  on;

%Obliczanie spektrum amplitudowego
spectrum = fft(ASK_signal);
spectrum_amp = abs(spectrum);

%Wykres spektrum aplitudowe
%subplot(5,1,2);
axes(handles.axes2)
plot(spectrum_amp);
xlabel('Spectral Index[Hz]');
ylabel('amplituda widma');
title('Widmo amplitudowe sygna�u');
grid on;

%ODBIORNIK
%snr = -10;
% bandwidth = 1;
% 
% received_max_level = max(abs(ASK_signal));
% noise_level = received_max_level / (10^(-snr/10))    %TODO Zapyta� na konsultacjach, p�ki co dzia�a
% m = size(ASK_signal);
% noise=randn(noise_level,m(2));

%DODANIE SZUMU
ASK_signal_with_noise = awgn(ASK_signal, snr);

%FILTRACJA
% windowSize = 5; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 4;
% ASK_signal_filtered = filter(b, a, ASK_signal_with_noise);

d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',0.15,'DesignMethod','butter');
% y = filtfilt(d1,x);
ASK_signal_filtered = filtfilt(d1, ASK_signal_with_noise);


% Wykres sygna�u z szumem
% subplot(5,1,3);
% plot(time,ASK_signal_with_noise,'LineWidth',1);
% xlabel('Time[s]');
% ylabel('Amplitude');
% title('Przebieg sygna�u z szumem w funkcji czasu');
% %axis([0 time(end) 1.5 1.5]);
% grid  on

% Wykres sygna�u po filtracji szumu
% subplot(5,1,4);
axes(handles.axes3)
plot(time,ASK_signal_filtered,'LineWidth',1);
xlabel('Time[s]');
ylabel('Amplitude');
title('Przebieg sygna�u po filtracji w funkcji czasu');
%axis([0 time(end) 1.5 1.5]);
grid  on

bits_number = length(ASK_signal_filtered)/ length(t);
ASK_signal_splitted = reshape(ASK_signal_filtered, [length(t), bits_number]);
ASK_max_values = max(ASK_signal_splitted);
ASK_signal_decoded = round(ASK_max_values);

%subplot(5,1,5);
axes(handles.axes4)
scatter([1:length(ASK_signal_decoded)], ASK_signal_decoded);
xlabel('Bit');
ylabel('Warto��');
title('Odebrany sygna� bitowy');
%axis([0 time(end) 1.5 1.5]);
grid  on

ASK_signal_chars = reshape(ASK_signal_decoded, 8, length(ASK_signal_decoded)/8);
ASK_signal_chars=ASK_signal_chars.';

errorArray = mod(sum(ASK_signal_chars,2),2);
% ASK_signal_chars(:,1) = 0;

% ASK_signal_chars
% ASK_signal_decoded
% X = num2str(ASK_signal_decoded);
% X = X(~isspace(X))
ASK_signal_chars
received_text = char(bin2dec(num2str(ASK_signal_chars))).';
received_text
sum(errorArray)

%display = sprintf('%s', received_text);
handles.text8.String=received_text

handles.text9.String=sum(errorArray)

