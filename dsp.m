% Exp 1: Signal Generation
clc; clear; close all;
t = -1:0.001:1; n = -5:5; m = -10:10; f = 5; A = 2;
T = {t,t,t,t,t,t,0:0.1:1,-5:0.001:5};
X = {double(t==0),double(t>=0),square(2*pi*f*t),sawtooth(2*pi*f*t),A*sawtooth(2*pi*f*t,0.25),sin(2*pi*f*t),0:2:20,sinc(-5:0.001:5)};
name = {'Unit impulse','Unit step','Square wave','Sawtooth wave','Triangular wave','Sinusoidal wave','Ramp signal','Sinc function'};
figure;
for i = 1:8
    subplot(4,2,i); plot(T{i},X{i}); title(name{i});
    xlabel('Time'); ylabel('Amplitude'); grid on;
end
N = {n,n,n,n,m,m,0:10,m};
Y = {double(n==0),double(n>=0),square(2*pi*f*n),sawtooth(pi*n/2),A*sawtooth(pi*m/4,0.5),sin(pi*m/4),0:2:20,sinc(m/2)};
figure;
for i = 1:8
    subplot(4,2,i); stem(N{i},Y{i},'filled'); title(name{i});
    xlabel('Time'); ylabel('Amplitude'); grid on;
end

% Exp 2: Discrete-Time Signal Operations
clc; clear; close all;
n = -4:4; x1 = 0:8; x2 = [2 3 4 5 6 7 0 0 1]; k = 2;
add = x1+x2; sub = x1-x2; mul = x1.*x2;
scale1 = 2*x1; scale2 = 0.5*x2;
shift1 = x1+2; shift2 = x2+2;
delay = [zeros(1,k) x1(1:end-k)];
advance = [x1(k+1:end) zeros(1,k)];
figure;
subplot(3,2,1); stem(n,x1,'filled'); hold on; stem(n,x2,'r'); stem(n,add,'g'); title('Addition');
subplot(3,2,2); stem(n,x1,'filled'); hold on; stem(n,x2,'r'); stem(n,sub,'g'); title('Subtraction');
subplot(3,2,3); stem(n,x1,'filled'); hold on; stem(n,x2,'r'); stem(n,mul,'g'); title('Multiplication');
subplot(3,2,4); stem(n,scale1,'filled'); hold on; stem(n,scale2,'r'); title('Amplitude scaling');
subplot(3,2,5); stem(n,shift1,'filled'); hold on; stem(n,shift2,'r'); title('Amplitude shifting');
subplot(3,2,6); stem(n,x1,'filled'); hold on; stem(n,delay,'r'); stem(n,advance,'g'); title('Time delay and advance');
for i = 1:6
    subplot(3,2,i); xlabel('n'); ylabel('Amplitude'); grid on;
end
subplot(3,2,1); legend('x1','x2','Addition');
subplot(3,2,2); legend('x1','x2','Subtraction');
subplot(3,2,3); legend('x1','x2','Multiplication');
subplot(3,2,4); legend('2*x1','0.5*x2');
subplot(3,2,5); legend('x1+2','x2+2');
subplot(3,2,6); legend('x1','Delayed x1','Advanced x1');
disp(add); disp(sub); disp(mul);
disp(scale1); disp(scale2); disp(shift1); disp(shift2);
disp(delay); disp(advance);

% Exp 3: Convolution and Cross-Correlation
clc; clear; close all;
x = [1 2 3 4]; h = [6 7 8];
y = conv(x,h);
r = conv(x,fliplr(h));
n = 0:length(y)-1;
disp('Linear convolution:'); disp(y);
disp('Cross-correlation:'); disp(r);
figure;
subplot(2,1,1); stem(n,y,'filled'); title('Linear convolution');
xlabel('n'); ylabel('Amplitude'); grid on;
subplot(2,1,2); stem(n,r,'filled'); title('Cross-correlation');
xlabel('n'); ylabel('Amplitude'); grid on;

% Exp 4: DFT and IDFT
clc; clear; close all;
x = [1 2 3 1];
N = length(x); n = 0:N-1; X = zeros(1,N);
for k = 0:N-1
    X(k+1) = sum(x .* exp(-1j*2*pi*k*n/N));
end
xr = zeros(1,N);
for k = 0:N-1
    xr(k+1) = sum(X .* exp(1j*2*pi*k*n/N))/N;
end
disp('DFT:'); disp(X);
disp('IDFT:'); disp(real(xr));
figure;
subplot(3,1,1); stem(n,x,'filled'); title('Input signal x[n]');
xlabel('n'); ylabel('Amplitude'); grid on;
subplot(3,1,2); stem(n,abs(X),'filled'); title('DFT magnitude |X[k]|');
xlabel('k'); ylabel('Magnitude'); grid on;
subplot(3,1,3); stem(n,angle(X),'filled'); title('Phase of DFT');
xlabel('k'); ylabel('Phase (rad)'); grid on;

% Exp 5: FFT (8-Point Radix-2 DIT-FFT)
clc; clear; close all;
x = [1 4 3 2 0 0 0 0]; N = 8;
X = x([1 5 3 7 2 6 4 8]);
for s = 1:3
    m = 2^s; w = exp(-1j*2*pi/m);
    for p = 1:m:N
        for j = 0:m/2-1
            a = X(p+j); b = X(p+j+m/2)*w^j;
            X(p+j) = a+b; X(p+j+m/2) = a-b;
        end
    end
end
disp('8-point DFT:'); disp(X);
figure;
subplot(2,1,1); stem(0:N-1,x,'filled'); title('Input sequence x[n]');
xlabel('n'); ylabel('Amplitude'); grid on;
subplot(2,1,2); stem(0:N-1,abs(X),'filled'); title('Magnitude of 8-point DFT');
xlabel('k'); ylabel('|X[k]|'); grid on;

% Exp 6: Convolution Theorem
clc; clear; close all;
x = [1 2 3 4]; h = [1 -1 1];
y = conv(x,h); N = length(y);
yf = ifft(fft(x,N).*fft(h,N));
disp('Time-domain convolution:'); disp(y);
disp('Frequency-domain result:'); disp(real(yf));
if max(abs(y-real(yf))) < 1e-10
    disp('Verification passed: Results are identical.');
else
    disp('Verification failed: Results differ.');
end
figure;
stem(0:N-1,y,'filled'); hold on;
plot(0:N-1,real(yf),'r--x');
title('Convolution theorem verification'); xlabel('n'); ylabel('Amplitude');
legend('Time-domain convolution','IFFT(FFT(x).*FFT(h))'); grid on;

% Exp 7: FIR Low-Pass Filter
clc; clear; close all;
N = 25; w = 0:0.01:pi;
win = {hamming(N+1),boxcar(N+1),bartlett(N+1),blackman(N+1),hanning(N+1)};
name = {'Hamming','Boxcar','Bartlett','Blackman','Hanning'};
figure; hold on;
for i = 1:5
    b = fir1(N,0.5,win{i});
    plot(w/pi,abs(freqz(b,1,w)));
end
xlabel('Normalized Frequency'); ylabel('|H(w)|');
title('FIR low-pass filter using window techniques');
legend(name); grid on;

% Exp 8: IIR Butterworth Low-Pass Filter
clc; clear; close all;
Fs = 10000; Wn = 3000/(Fs/2); order = [4 7 9];
figure;
subplot(2,1,1); hold on; grid on; title('Butterworth magnitude response');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); ylim([-0.5 1.5]);
subplot(2,1,2); hold on; grid on; title('Magnitude response in dB');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); ylim([-80 20]);
for i = 1:3
    [b,a] = butter(order(i),Wn,'low');
    [H,f] = freqz(b,a,100,Fs);
    subplot(2,1,1); plot(f,abs(H));
    subplot(2,1,2); plot(f,20*log10(abs(H)));
end
subplot(2,1,1); legend('Order 4','Order 7','Order 9');
subplot(2,1,2); legend('Order 4','Order 7','Order 9');

% Exp 9: ECG Signal Generation
clc; clear; close all;
t = 0:0.005:1;
P = 0.2*sin(2*pi*5*(t-0.05)).*(t>=0.05 & t<=0.15);
QRS = -1.2*sawtooth(2*pi*20*(t-0.3)).*(t>=0.3 & t<=0.35);
T = 0.5*sin(2*pi*5*(t-0.6)).*(t>=0.6 & t<=0.7);
ecg = P + QRS + T;
figure; plot(t,ecg); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('One ECG Cycle');
y = repmat(ecg,1,5); ty = linspace(0,5,length(y));
figure; plot(ty,y); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('ECG Signal');

% Exp 10: Musical Tone Generation
clc; clear; close all;
Fs = 44100; notes = [240 240 270 270 300 300 320 320];
t = 0:1/Fs:0.5; signal = [];
for f = notes
    signal = [signal 0.5*sin(2*pi*f*t)];
end
sound(signal,Fs);
total_time = (0:length(signal)-1)/Fs;
figure; plot(total_time,signal,'LineWidth',1.1);
xlabel('Time (seconds)'); ylabel('Amplitude'); grid on;
title('Musical Tone Sequence: Saa Saa Re Re Ga Ga Ma Ma');
k = 1:round(Fs*0.01);
figure; plot(total_time(k),signal(k),'LineWidth',1.2);
xlabel('Time'); ylabel('Amplitude'); grid on;
title('Zoomed View of Single Tone Generation');
