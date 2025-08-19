# EEG Oddball Analysis for N-word

## Overview
This MATLAB script analyzes EEG data from an oddball experiment involving 5 African American male participants (PT1NN, PT2YY, PT3YY, PT4YY, PT7YN) to investigate neural responses to the N-word. It includes preprocessing, artifact rejection, ICA, epoching, normalization, peak extraction, and computes metrics like Frontal Alpha Asymmetry (FAA), Beta/Alpha Ratio (BAR), Theta/Beta Ratio (TBR), and High Beta Power (HB). Outputs include tables, plots, topoplots for hypotheses H1-H4, and PCA.

## Requirements
- MATLAB (tested on R2023a+)
- EEGLAB toolbox (v2025.0.0): Download from [EEGLAB website](https://sccn.ucsd.edu/eeglab/downloadtoolbox.php) and place the extracted `eeglab2025.0.0` folder in the project directory.
- Input files: `<participant_id>.csv`, `<participant_id>_Interval Marker.csv` for each participant
- Channels: AF3, AF4, F3, F4, FC5, FC6, P8, F8

## Usage
1. Place input CSV files in the working directory.
2. Download and extract the EEGLAB v2025.0.0 folder from the EEGLAB website and place it in the project directory.
3. Run the script: `run('Ner_word_EEG_analysis.m')`
4. Outputs: CSV tables, PNG plots, MAT files (e.g., PeakAnalysisResults.mat)

## Sections
- Data Extraction: Loads CSV data/markers.
- Preprocessing: Bandpass filter (0.8-40 Hz).
- Artifact Rejection: Joint probability, amplitude threshold (±100 µV).
- ICA: Multi-channel ICA, manual artifact removal (commented).
- Epoching: -200 to 1000 ms around events ('W', 'T1', 'T2').
- Normalization: Log-transform, z-score.
- Peak Extraction: Peaks in time windows (P200, N200, etc.) and bands (Theta, Alpha, Beta, HighBeta).
- Metrics Calculation:
  - FAA: F4 - F3 alpha power.
  - BAR: Beta/Alpha ratio.
  - TBR: Theta/Beta ratio.
  - HB: High Beta power.
- Hypotheses Analysis (H1-H4): Tables, plots, stats (Wilcoxon, Spearman's, Cohen's d).
- PCA: On metrics for T1/T2 conditions.
- Topoplots: HB and FAA for W, T1, T2, T1T2.

## Outputs
- CSVs: Segment tables (FAA, BAR, etc.), stats, PCA.
- PNGs: Peak/Z-score graphs, topoplots, bar charts per participant/group.

## Notes
- Assumes fs=256 Hz.
- Manual ICA removal commented; adjust as needed.
- Y* group: PT2YY, PT3YY, PT4YY, PT7YN.

## Credits
This work belongs to Regina Brown. I assisted with the EEG analysis. Her dissertation defense: "Language, Emotion, and Cognitive Congruence: Does Appropriated Racism Detoxify Ner/Na for African American Males?" on July 22, 2025, 11:00 AM CT, Zoom or Saint Mary's Minneapolis Campus. RSVP: https://calendly.com/regina2-0/regina-brown-dissertation-defense-rsvp.

## License
[MIT License](#license)

## License
MIT License

Copyright (c) 2025 Jesklose

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.