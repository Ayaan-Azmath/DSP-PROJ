# Harmonic-Percussive Source Separation Using STFT and Median Filtering

MATLAB implementation of harmonic-percussive source separation (HPSS) for music signals, based on median filtering of the STFT spectrogram (Fitzgerald, 2010). The project studies how mask type, median-filter length, window length and additive noise affect separation.

**Author:** Mohammed Ayaan Azmath (Roll No. 1601-24-735-050)
**Supervisor:** Dr. Bhasker Dappuri
**Institution:** Chaitanya Bharathi Institute of Technology, Hyderabad, India
**Course project:** Digital Signal Processing

## Method

1. Zero-pad the signal and compute the STFT, then the magnitude spectrogram `S`.
2. Median-filter `S` along time (harmonic-enhanced `H`) and along frequency (percussive-enhanced `P`).
3. Build binary or soft masks from `H` and `P`.
4. Apply the masks to the complex STFT and reconstruct both components with the inverse STFT.
5. Remove the padding. The padding avoids large spikes at the first and last samples of the inverse transform.

## Files

| File | Purpose |
|---|---|
| `hpss_main.m` | Runs HPSS on one audio clip. Plots spectrograms, masks and waveforms, saves PNG figures and the separated WAV files, and prints energy shares and reconstruction error. |
| `hpss_experiments.m` | Synthetic benchmark with known references (SDR for mask type, filter length, window length and noise) and a real-clip energy table. Prints results as LaTeX table rows. |
| `figures/` | Output figures (spectrograms, masks, waveforms). |

## Requirements

- MATLAB R2019a or later
- Signal Processing Toolbox (`stft`, `istft`)

## How to run

1. Put your audio clip in the same folder as the scripts. Audio files are not included in this repository.
2. In `hpss_main.m` and `hpss_experiments.m`, set `audioFile` to your file name (for example `music.wav`).
3. Set MATLAB's Current Folder to this folder.
4. Run `hpss_main.m` for the separation and figures, or `hpss_experiments.m` for the benchmark tables.

The synthetic benchmark in `hpss_experiments.m` needs no audio file. The random seed is fixed at 2026.

## Default parameters

| Parameter | Value |
|---|---|
| Window | Hann, 1024 samples |
| Hop size | 256 samples |
| Median filter length (time and frequency) | 17 bins |
| Mask | soft |

## Results (summary)

Real excerpt, 10 s piano clip at 44.1 kHz, soft masks: 89.2% of the energy in the harmonic component, 1.5% in the percussive component, reconstruction error about 1e-16.

Synthetic benchmark (SDR in dB, harmonic / percussive):

| Setting | SDR_H | SDR_P |
|---|---|---|
| Mixture used as estimate | 13.14 | -13.14 |
| Binary mask (L=17, N=1024) | 21.24 | 8.10 |
| Soft mask (L=17, N=1024) | 15.09 | 1.94 |
| Soft, L=51 | 20.89 | 7.74 |
| Soft, N=2048 | 21.11 | 7.97 |

## Limitations

- One real clip without reference stems, so real-clip results use energy shares and visual inspection only.
- The synthetic benchmark is simple (stationary chords and drum-like bursts).
- SDR is the simple form, not the full BSS Eval measure.
- The sweeps over filter length and window length did not reach an optimum.
- No listening tests and no comparison with other HPSS or learned methods.

## Reference

D. Fitzgerald, "Harmonic/percussive separation using median filtering," in *Proc. 13th Int. Conf. Digital Audio Effects (DAFx)*, Graz, Austria, 2010.
