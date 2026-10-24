%% HPSS experiments: synthetic benchmark (SDR) + real-clip energy table
% Prints rows in LaTeX table format so you can paste them into main_v2.tex.
% Requires: Signal Processing Toolbox (stft/istft), MATLAB R2019a or later.
% Optional: put your real clip in this folder as music.wav for Table III.

clc; clear; close all; rng(2026);

%% Settings
fs = 22050;  dur = 10;
N0 = 1024;   H0 = 256;   L0 = 17;

[x, xh0, xp0] = makeSynthetic(fs, dur);

%% Baseline
fprintf('=== Synthetic benchmark (SDR in dB) ===\n');
fprintf('Mixture used as estimate & %.2f & %.2f \\\\\n', sdr(xh0, x), sdr(xp0, x));

%% Exp 1: mask type
fprintf('\n-- Mask type (N=%d, H=%d, L=%d) --\n', N0, H0, L0);
for mt = {'binary', 'soft'}
    [sH, sP] = evalCfg(x, xh0, xp0, fs, N0, H0, L0, L0, mt{1});
    fprintf('%s mask & %.2f & %.2f \\\\\n', mt{1}, sH, sP);
end

%% Exp 2: median filter length
fprintf('\n-- Filter length L (soft, N=%d, H=%d) --\n', N0, H0);
for L = [5 9 17 31 51]
    [sH, sP] = evalCfg(x, xh0, xp0, fs, N0, H0, L, L, 'soft');
    fprintf('$L=%d$ & %.2f & %.2f \\\\\n', L, sH, sP);
end

%% Exp 3: window length (hop = N/4)
fprintf('\n-- Window length N (soft, L=%d, H=N/4) --\n', L0);
for N = [512 1024 2048]
    [sH, sP] = evalCfg(x, xh0, xp0, fs, N, N/4, L0, L0, 'soft');
    fprintf('$N=%d$ & %.2f & %.2f \\\\\n', N, sH, sP);
end

%% Exp 4: additive noise (references stay clean; noise counts as distortion)
fprintf('\n-- Additive white noise (soft, N=%d, H=%d, L=%d) --\n', N0, H0, L0);
for snr = [Inf 20 10 0]
    if isinf(snr)
        xn = x;
    else
        Pn = mean(x.^2) / 10^(snr/10);
        xn = x + sqrt(Pn) * randn(size(x));
    end
    [xh, xp] = hpss(xn, fs, N0, H0, L0, L0, 'soft');
    fprintf('SNR $=%g$ dB & %.2f & %.2f \\\\\n', snr, sdr(xh0, xh), sdr(xp0, xp));
end

%% Real clip: energy share and reconstruction error (Table III)
fprintf('\n=== Real clip (energy %%, reconstruction error) ===\n');
audioFile = 'music.wav';   % <-- change this to your clip's exact file name (or full path)
if exist(audioFile, 'file')
    [xr, fsr] = audioread(audioFile);
    xr = mean(xr, 2);
    xr = xr(1:min(length(xr), round(10*fsr)));
    xr = xr / max(abs(xr));
    cfgs = {'soft', 17; 'binary', 17; 'soft', 9; 'soft', 31};
    Ex = sum(xr.^2);
    for k = 1:size(cfgs, 1)
        Lk = cfgs{k, 2};
        [xh, xp] = hpss(xr, fsr, 1024, 256, Lk, Lk, cfgs{k, 1});
        fprintf('%s, $L=%d$ & %.1f & %.1f & %.2e \\\\\n', cfgs{k, 1}, Lk, ...
            100*sum(xh.^2)/Ex, 100*sum(xp.^2)/Ex, sum((xr - xh - xp).^2)/Ex);
    end
else
    fprintf('"%s" not found. MATLAB is looking in: %s\n', audioFile, pwd);
    fprintf('Copy your clip into that folder, or set audioFile to its full path.\n');
end

%% ===================== Local functions =====================

function [sH, sP] = evalCfg(x, xh0, xp0, fs, N, H, Lh, Lp, mt)
    [xh, xp] = hpss(x, fs, N, H, Lh, Lp, mt);
    sH = sdr(xh0, xh);
    sP = sdr(xp0, xp);
end

function s = sdr(ref, est)
    ref = ref(:); est = est(:);
    s = 10*log10(sum(ref.^2) / sum((ref - est).^2));
end

function [xh, xp] = hpss(x, fs, N, H, Lh, Lp, maskType)
    x = x(:);  L = length(x);
    win = hann(N, 'periodic');  ov = N - H;
    xpad = [zeros(N,1); x; zeros(2*N,1)];          % padding avoids istft edge spikes
    X = stft(xpad, fs, 'Window', win, 'OverlapLength', ov, ...
             'FFTLength', N, 'FrequencyRange', 'onesided');
    S  = abs(X);
    Hh = movmedian(S, Lh, 2);                      % along time
    Pp = movmedian(S, Lp, 1);                      % along frequency
    switch lower(maskType)
        case 'binary'
            Mh = double(Hh >= Pp);  Mp = 1 - Mh;
        otherwise
            e = 1e-10;
            Mh = Hh.^2 ./ (Hh.^2 + Pp.^2 + e);
            Mp = Pp.^2 ./ (Hh.^2 + Pp.^2 + e);
    end
    xh = istft(Mh .* X, fs, 'Window', win, 'OverlapLength', ov, ...
               'FFTLength', N, 'ConjugateSymmetric', true, 'FrequencyRange', 'onesided');
    xp = istft(Mp .* X, fs, 'Window', win, 'OverlapLength', ov, ...
               'FFTLength', N, 'ConjugateSymmetric', true, 'FrequencyRange', 'onesided');
    xh = fixLen(real(xh(N+1:end)), L);
    xp = fixLen(real(xp(N+1:end)), L);
end

function y = fixLen(y, L)
    y = y(:);
    if length(y) >= L
        y = y(1:L);
    else
        y = [y; zeros(L - length(y), 1)];
    end
end

function [x, xh, xp] = makeSynthetic(fs, dur)
    n = round(dur*fs);
    t = (0:n-1)'/fs;

    % Harmonic part: four chords, five harmonics per note
    chords = [220 277.18 329.63; 196 246.94 293.66; 174.61 220 261.63; 196 246.94 293.66];
    segLen = round(n / size(chords, 1));
    xh = zeros(n, 1);
    for c = 1:size(chords, 1)
        idx = (c-1)*segLen + (1:segLen);
        idx = idx(idx <= n);
        tt  = t(idx) - t(idx(1));
        env = (1 - exp(-tt/0.05)) .* exp(-0.3*tt);
        seg = zeros(size(tt));
        for f0 = chords(c, :)
            for h = 1:5
                seg = seg + (1/h) * sin(2*pi*h*f0*tt);
            end
        end
        xh(idx) = seg .* env;
    end
    xh = 0.6 * xh / max(abs(xh));

    % Percussive part: kick / snare / hi-hat on a 0.25 s grid
    tk = (0:round(0.15*fs)-1)'/fs;   kick  = sin(2*pi*60*tk) .* exp(-tk/0.04);
    ts = (0:round(0.20*fs)-1)'/fs;   snare = randn(size(ts)) .* exp(-ts/0.04);
    th = (0:round(0.05*fs)-1)'/fs;   hat   = [0; diff(randn(length(th)-1+1, 1))];
    hat = hat(1:length(th)) .* exp(-th/0.01);
    xp = zeros(n, 1);
    for k = 0:round(dur/0.25)-1
        i0 = round(k*0.25*fs) + 1;
        xp = addBurst(xp, i0, 0.3*hat);
        if mod(k, 4) == 0, xp = addBurst(xp, i0, 0.9*kick);  end
        if mod(k, 4) == 2, xp = addBurst(xp, i0, 0.5*snare); end
    end
    xp = 0.5 * xp / max(abs(xp));

    x  = xh + xp;
    sc = max(abs(x));
    x = x / sc;  xh = xh / sc;  xp = xp / sc;      % common scale keeps x = xh + xp
end

function y = addBurst(y, i0, b)
    m = min(length(b), length(y) - i0 + 1);
    if m > 0
        y(i0:i0+m-1) = y(i0:i0+m-1) + b(1:m);
    end
end
