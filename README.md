# Harmonic-Percussive Source Separation Using STFT and Median Filtering

MATLAB implementation of harmonic-percussive source separation (HPSS) for music signals, based on median filtering of the STFT spectrogram (Fitzgerald, 2010). The project studies how mask type, median-filter length, window length and additive noise affect separation.

- **Author:** Mohammed Ayaan Azmath (Roll No. 1601-24-735-050)
- **Supervisor:** Dr. Bhasker Dappuri
- **Institution:** Chaitanya Bharathi Institute of Technology, Hyderabad, India
- **Course project:** Digital Signal Processing

## Repository structure

```
.
├── code/
│   ├── hpss_main.m          # HPSS on one clip: plots, figures, separated audio
│   └── hpss_experiments.m   # synthetic SDR benchmark + real-clip energy table
├── paper/
│   ├── main.tex             # IEEE conference paper (LaTeX)
│   └── README.md
├── results/
│   ├── results_summary.md   # tables of all reported results
│   └── (figures written here by hpss_main.m)
├── requirements.txt
├── .gitignore
└── README.md
```

## Method

1. Zero-pad the signal and compute the STFT, then the magnitude spectrogram `S`.
2. Median-filter `S` along time (harmonic-enhanced `H`) and along frequency (percussive-enhanced `P`).
3. Build binary or soft masks from `H` and `P`.
4. Apply the masks to the complex STFT and reconstruct both components with the inverse STFT.
5. Remove the padding. The padding avoids large spikes at the first and last samples of the inverse transform.

## Requirements

See `requirements.txt`: MATLAB R2019a or later with the Signal Processing Toolbox.

## How to run

1. Put your audio clip in the `code/` folder. Audio files are not included in this repository.
2. In `hpss_main.m` and `hpss_experiments.m`, set `audioFile` to your file name (for example `music.wav`).
3. Set MATLAB's Current Folder to `code/`.
4. Run `hpss_main.m` for the separation. Figures and separated audio are saved to `results/`.
5. Run `hpss_experiments.m` for the benchmark tables. It prints results as LaTeX table rows. The synthetic benchmark needs no audio file, and the random seed is fixed at 2026.

## Default parameters

| Parameter | Value |
|---|---|
| Window | Hann, 1024 samples |
| Hop size | 256 samples |
| Median filter length (time and frequency) | 17 bins |
| Mask | soft |

## Results (summary)

Real excerpt (10 s piano clip, 44.1 kHz, soft masks): 89.2% of the energy in the harmonic component, 1.5% in the percussive component, reconstruction error about 1e-16.

Synthetic benchmark (SDR in dB, harmonic / percussive): using the mixture as the estimate gives 13.14 / -13.14. Soft masks give 15.09 / 1.94, and binary masks give 21.24 / 8.10. SDR increased with filter length and window length over the tested ranges. All tables are in `results/results_summary.md`.

## Limitations

- One real clip without reference stems, so real-clip results use energy shares and visual inspection only.
- The synthetic benchmark is simple (stationary chords and drum-like bursts).
- SDR is the simple form, not the full BSS Eval measure.
- The sweeps over filter length and window length did not reach an optimum.
- No listening tests and no comparison with other HPSS or learned methods.

## Reference

D. Fitzgerald, "Harmonic/percussive separation using median filtering," in *Proc. 13th Int. Conf. Digital Audio Effects (DAFx)*, Graz, Austria, 2010.
