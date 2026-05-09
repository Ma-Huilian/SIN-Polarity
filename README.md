# SIN-Polarity

A Novel Algorithm for Automatic Identification of P-wave First-motion Polarity Based on a Sine-wave Template.

## Overview

This project implements a workflow for:
1. Fetching seismic waveform data from IRIS
2. Picking P-wave arrival times (PhaseNet or STA/LTA)
3. Determining first motion polarity by calculating the cross-correlation values between P-wave and sine wave template

## Examples

### Example 1: Oklahoma Earthquake (2014-10-10, mb4.5)

### Example 2: Myanmar Earthquake (2025-03-28, mw7.7)

## Workflow

```
SIN_Polarity_Step1_Fetch_data.py                # Step 1: Download data from IRIS
        ↓
SIN_Polarity_Step2_First_arrival_picking.py     # Step 2: Pick P-wave first arrival by PhaseNet or STA/LTA
        ↓
SIN_Polarity_Step3_Polarity_identification.py   # Step 3: Determine P-wave first-motion polarity
```

## Quick Start

### Step 1: Download Waveform Data

Run `SIN_Polarity_Step1_Fetch_data.py` to download seismic waveform data from IRIS.

Key parameters to configure:
- `origin_time`: Event origin time (UTC)
- `lat`, `lon`: Epicenter coordinates
- `minDisInDeg`, `maxDisInDeg`: Station distance range (in degrees)

### Step 2: Pick P-wave First Arrival

Run `SIN_Polarity_Step2_First_arrival_picking.py` to pick P-wave arrival times using PhaseNet (example 1) or STA/LTA (example 2).

Before running:
1. Organize three-component data into separate folders: `data1/` (E), `data2/` (N), `dataZ/` (Z)
2. Set event location parameters: `ev_lat`, `ev_lon`, `ev_depth_km`
3. Prepare instrument response files (XML format)

Output: `templates.csv` containing P-wave arrival times for each station.

### Step 3: Determine P-wave First-Motion Polarity

Run `SIN_Polarity_Step3_Polarity_identification.py` to determine P-wave first-motion polarity on the vertical component.

Input:
- `templates.csv`: P-wave arrival times
- `dataENZ/`: Three-component waveform data (MSEED format)

Output: `StationinformationPhasenet.csv` / `StationinformationSTALTA.csv`: polarity results (1='U', -1='D').

### Quick Path: Directly Use Step 3

If you already have P-wave first arrival times from public datasets or other picking methods, you can directly use `SIN_Polarity_Step3_Polarity_identification.py` to determine P-wave first-motion polarity.

Prepare your arrival time data in `templates.csv` format with the following columns:
- `network`: Network code
- `station`: Station name
- `onsetP`: P-wave arrival time (seconds from trace start)

Example:
```
network,station,onsetP
GS,OK025,41.54
NX,STN02,41.9
```

## Requirements

- Python 3.x
- obspy
- numpy
- pandas
- scipy
- tensorflow (for PhaseNet)
- matplotlib
- numba (for focal mechanism inversion)

## Contact

Any questions or advices? Please contact:
- huilianma@stu.cdut.edu.cn
- jianlongyuan@cdut.edu.cn

