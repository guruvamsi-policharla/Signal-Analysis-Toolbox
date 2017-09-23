% Plot wavelet transform

sig=LWBF(1,1:4:end); % Resampled for speed
fs=10;
fmin=0.005;
fmax=2;
f0=1;

time=linspace(0,length(sig)/fs,length(sig));

[WT,freq]=wt(sig,fs,'fmin',fmin,'fmax',fmax,'CutEdges','on','f0',f0);

figure
surf(time,freq,abs(WT))
shading interp
view(0,90)
set(gca,'yscale','log')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
axis tight