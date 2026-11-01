# SH Project: PMT pulse analysis

## Folders

| Folder | Contents |
|---|---|
| `code/` | Current analysis macros. Keep this in sync with the lab PC. |
| `raw_data/` | Digitiser `.dat` files (CAEN WaveDump binary: 24-byte header + 1024 float32 samples per event). Never write here. |
| `results/` | Macro outputs, one `<run>_results/` folder per raw file, plus `dark_count/`. Can be regenerated. |
| `figures/` | Plots picked for lab notes and the report. |
| `archive/` | Old Python scripts and `read_wave.cpp`, kept for reference only. |
| `literature/` | Papers. |

## Running

`pmt_events.C` saves to `<data folder>_results/<data file>/` next to the macro, with files named
`<data folder>_<input name>_pulse_*`. Example on the lab PC:

```
root -l 'pmt_events.C("D:/data/23_10_2026/wave_0.dat")'
#  -> <macro folder>/23_10_2026_results/wave_0/23_10_2026_wave_0_pulse_events.root, _pulse_stats, _pulse_times, _pulse_waveforms
#     (each plot as .png, plus .pdf when SAVE_PDF = true)
```

On the laptop, pass `../results` as the third argument so outputs go into `results/` rather than `code/`:

```
cd code
root -l 'pmt_events.C("../raw_data/wave_0.dat", -1, "../results")'   # -> results/raw_data_results/wave_0/
root -l dark_count.C                                                 # -> results/dark_count/
```

## Pulse viewer

`pulse_viewer.C` shows the pulses one at a time. It `#include`s `pmt_events.C` (keep both in the same
folder) and uses its settings and `analyse()`, so it finds exactly the same pulses: same `DATA_FILE`,
`THRESH_SIGMA` and time gate. Run it in a fresh ROOT session, not after `.L pmt_events.C`.

```
cd code
root -l pulse_viewer.C                              # pulses in the time gate (same as the tree)
root -l 'pulse_viewer.C("D:/data/23_10_2026/wave_0.dat")'   # another data file
root -l 'pulse_viewer.C(DATA_FILE, 3.0)'            # only pulses taller than 3 mV
root -l 'pulse_viewer.C(DATA_FILE, 0, 488)'         # only pulses peaking in bin 488 (97.6 ns)
root -l 'pulse_viewer.C(DATA_FILE, 0, -1, true)'    # tallest point anywhere in the record, no time gate
```

Arguments: data file, minimum height (mV), peak bin (-1 = any), whole record instead of the time gate (true/false).
Don't use `-b` or `-q`: the viewer needs its windows and the ROOT prompt.

**Moving between pulses**

| Action | "Pulses" button window | Keys (click inside the plot first) | ROOT prompt |
|---|---|---|---|
| next / previous | `next >` / `< previous` | `n` / `p` or arrows | `nextPulse()` / `prevPulse()` |
| +100 / -100 | `+100 >>` / `<< -100` | `N` / `P` | `jumpPulses(100)` / `jumpPulses(-100)` |
| any jump | | | `jumpPulses(n)` |
| pulse number k of the selection | | | `showPulse(k)` |
| event number from the file / tree | | | `showEvent(95056)` |
| save the plot as .png | `save png` | `s` | `savePulse()` |

Saved plots go to `<data folder>_results/<data file>/pulses/event_<n>.png` next to the macro. Nothing
else is written.

**What the plots show**

- Top: the whole 205 ns waveform. Yellow = time gate (85-115 ns with the current settings, `GATE_CENTER_NS` +/- `GATE_WIDTH_NS`/2,
  the only place the pulse is looked for); blue = baseline (2-42 ns average);
  red = threshold (baseline - `THRESH_SIGMA` x rms); red dot = peak. Top right: empty events of the
  whole file (no pulse anywhere in 0-204.8 ns, time gate not applied).
- Bottom: zoom from 5 ns before to 15 ns after the peak, every 0.2 ns sample marked.
  Green = integration window, i.e. the samples summed for the charge. Title: height (mV and sigma),
  rms, peak time, charge Q and integration window length.

**Terminal output**: data file processed, threshold in sigma, number of events scanned and pulses selected,
empty events, then one line per pulse shown (event, height, peak time and bin, Q).

**Settings** at the top of `pulse_viewer.C`:

| Setting | Default | Meaning |
|---|---|---|
| `VIEW_ALL` | `false` | `true`: tallest point anywhere in the record (0-204.8 ns), not only the time gate |
| `MIN_HEIGHT_MV` | `0` | only pulses at least this tall |
| `ONLY_PEAK_BIN` | `-1` | only pulses peaking in this bin (-1 = any) |
| `SORT_BY_HEIGHT` | `false` | `true`: tallest pulse first |
| `ZOOM_BEFORE_NS`, `ZOOM_AFTER_NS` | `5`, `15` | bottom plot range around the peak |
| `FIXED_Y` | `true` | same voltage range for every pulse, relative to its baseline (`false` = fit each waveform) |
| `Y_BELOW_MV`, `Y_ABOVE_MV` | `8`, `5` | fixed range: baseline - 8 mV to baseline + 5 mV (anything outside is cut off) |

`MIN_HEIGHT_MV`, `ONLY_PEAK_BIN` and `VIEW_ALL` are the defaults of the arguments above; an argument overrides them.

## Run log

| Raw file | Date | HV (V) | Setup (dark / LED / source) | Counting time | Notes |
|---|---|---|---|---|---|
| `wave_0.dat` | 2026-10-01 | | | | ~1.05M events |
| `wave_0__16_10_2026.dat` | 2026-10-16 | | | | ~347k events |
| dark count scan (-2 to -20 mV) | | | dark | | -4 mV read twice (26, 37); 26 used |
