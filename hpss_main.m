%% Harmonic-Percussive Source Separation using STFT and Median Filtering
% Requires: Signal Processing Toolbox (stft/istft), MATLAB R2019a or later.
% Usage: put a music file (e.g., music.wav) in the same folder and set
%        audioFile below. If the file is not found, a synthetic test
%        signal is generated so the script still runs.

clc; clear; close all;

%% 1. PARAMETERS
audioFile = 'music.wav';   % your music clip
maxDur    = 10;            % seconds to process
N         = 1024;          % window / FFT length
H         = 256;           % hop size
Lh        = 17;            % median filter length along time (harmonic)
Lp        = 17;            % median filter length along frequency (percussive)
maskType  = 'soft';        % 'soft' or 'binary'

%% 2. LOAD AUDIO (or generate synthetic signal)
synthetic = false;
if exist(audioFile, 'file')
    [x, fs] = audioread(audioFile);
    x = mean(x, 2);                          % convert to mono
    x = x(1:min(length(x), round(maxDur*fs)));
else
    synthetic = true;
    fs = 22050;
    t  = (0:round(maxDur*fs)-1)'/fs;
    % Harmonic part: sustained chords
    xh_true = 0.3*sin(2*pi*220*t) + 0.2*sin(2*pi*330*t) + ...
              0.15*sin(2*pi*440*t) + 0.1*sin(2*pi*660*t);
    xh_true = xh_true .* (0.7 + 0.3*sin(2*pi*0.5*t));
    % Percussive part: short noise bursts every 0.5 s
    rng(1);
    xp_true = zeros(size(t));
    burstLen = round(0.03*fs);
    env = exp(-(0:burstLen-1)'/(0.006*fs));
    for k = 1:round(maxDur/0.5)
        idx = round((k-1)*0.5*fs) + (1:burstLen);
        idx = idx(idx <= length(xp_true));
        xp_true(idx) = xp_true(idx) + 0.8*randn(length(idx),1).*env(1:length(idx));
    end
    x = xh_true + xp_true;
    fprintf('Audio file not found. Using a synthetic test signal.\n');
end
x = x / max(abs(x));                         % normalize
fprintf('Sampling rate: %d Hz, duration: %.2f s\n', fs, length(x)/fs);

%% 3. STFT
win     = hann(N, 'periodic');
overlap = N - H;
xpad = [zeros(N,1); x(:); zeros(2*N,1)];     % padding avoids istft edge artifacts
[X, f, tt] = stft(xpad, fs, 'Window', win, 'OverlapLength', overlap, ...
                  'FFTLength', N, 'FrequencyRange', 'onesided');
S = abs(X);                                  % magnitude spectrogram

%% 4. MEDIAN FILTERING
Hh = movmedian(S, Lh, 2);                    % along time  -> harmonic enhanced
Pp = movmedian(S, Lp, 1);                    % along freq  -> percussive enhanced

%% 5. MASKS
switch lower(maskType)
    case 'binary'
        Mh = double(Hh >= Pp);
        Mp = 1 - Mh;
    case 'soft'
        epsv = 1e-10;
        Mh = Hh.^2 ./ (Hh.^2 + Pp.^2 + epsv);
        Mp = Pp.^2 ./ (Hh.^2 + Pp.^2 + epsv);
end

%% 6. APPLY MASKS AND RECONSTRUCT (inverse STFT)
Xh = Mh .* X;
Xp = Mp .* X;

xh = istft(Xh, fs, 'Window', win, 'OverlapLength', overlap, ...
           'FFTLength', N, 'ConjugateSymmetric', true, 'FrequencyRange', 'onesided');
xp = istft(Xp, fs, 'Window', win, 'OverlapLength', overlap, ...
           'FFTLength', N, 'ConjugateSymmetric', true, 'FrequencyRange', 'onesided');
xh = real(xh); xp = real(xp);
xh = xh(N+1:end);                            % remove the padding added before STFT
xp = xp(N+1:end);

% Match lengths to the original
L  = length(x);
xh = fixLen(xh, L);
xp = fixLen(xp, L);

%% 7. EVALUATION
Ex  = sum(x.^2);
Eh  = sum(xh.^2);
Ep  = sum(xp.^2);
recErr = sum((x - (xh + xp)).^2) / Ex;
fprintf('\n--- Results (%s mask) ---\n', maskType);
fprintf('Harmonic energy   : %.1f %%\n', 100*Eh/Ex);
fprintf('Percussive energy : %.1f %%\n', 100*Ep/Ex);
fprintf('Reconstruction error (x - (xh+xp)): %.2e\n', recErr);

if synthetic
    % Ground truth is known only for the synthetic signal
    sc = max(abs(xh_true + xp_true));
    xh_ref = xh_true / sc;  xp_ref = xp_true / sc;
    sdr = @(ref, est) 10*log10(sum(ref.^2) / sum((ref - est).^2));
    fprintf('SDR harmonic   : %.2f dB\n', sdr(xh_ref, xh));
    fprintf('SDR percussive : %.2f dB\n', sdr(xp_ref, xp));
end

%% 8. PLOTS
dB = @(A) 20*log10(A + 1e-6);
fk = f/1000;  tsec = tt - N/fs;              % shift time axis back (padding removed)

fig1 = figure('Name','Spectrograms','Position',[100 100 900 750]);
subplot(3,1,1); imagesc(tsec, fk, dB(S)); axis xy; colorbar;
title('Original Spectrogram'); ylabel('Frequency (kHz)');
subplot(3,1,2); imagesc(tsec, fk, dB(abs(Xh))); axis xy; colorbar;
title('Harmonic Component'); ylabel('Frequency (kHz)');
subplot(3,1,3); imagesc(tsec, fk, dB(abs(Xp))); axis xy; colorbar;
title('Percussive Component'); ylabel('Frequency (kHz)'); xlabel('Time (s)');
for k = 1:3, subplot(3,1,k); clim([-80 max(dB(S(:)))]); xlim([0 L/fs]); end  % use caxis on older MATLAB

fig2 = figure('Name','Masks','Position',[150 150 900 500]);
subplot(2,1,1); imagesc(tsec, fk, Mh); axis xy; colorbar;
title('Harmonic Mask'); ylabel('Frequency (kHz)');
subplot(2,1,2); imagesc(tsec, fk, Mp); axis xy; colorbar;
title('Percussive Mask'); ylabel('Frequency (kHz)'); xlabel('Time (s)');

fig3 = figure('Name','Waveforms','Position',[200 200 900 600]);
tw = (0:L-1)'/fs;
subplot(3,1,1); plot(tw, x);  title('Original');   ylabel('Amplitude'); grid on;
subplot(3,1,2); plot(tw, xh); title('Harmonic');   ylabel('Amplitude'); grid on;
subplot(3,1,3); plot(tw, xp); title('Percussive'); ylabel('Amplitude'); xlabel('Time (s)'); grid on;

%% 9. SAVE OUTPUTS (images for the report, audio for the demo)
outDir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
exportgraphics(fig1, fullfile(outDir, 'spectrograms.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'masks.png'),        'Resolution', 300);
exportgraphics(fig3, fullfile(outDir, 'waveforms.png'),    'Resolution', 300);
audiowrite(fullfile(outDir, 'original.wav'),   x,  fs);
audiowrite(fullfile(outDir, 'harmonic.wav'),   xh/max(abs(xh)), fs);
audiowrite(fullfile(outDir, 'percussive.wav'), xp/max(abs(xp)), fs);
fprintf('\nSaved figures and wav files to: %s\n', outDir);

%% 10. PLAYBACK (uncomment to listen)
% sound(x,  fs); pause(length(x)/fs + 1);
% sound(xh, fs); pause(length(x)/fs + 1);
% sound(xp, fs);

%% Helper function
function y = fixLen(y, L)
    y = y(:);
    if length(y) >= L
        y = y(1:L);
    else
        y = [y; zeros(L - length(y), 1)];
    end
end
