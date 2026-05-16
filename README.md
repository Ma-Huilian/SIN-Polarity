# SIN-Polarity

A Novel Algorithm for Automatic Identification of P-wave First-motion Polarity Based on a Sine-wave Template.

## Overview

This project implements a workflow for:
1. Fetching seismic waveform data from IRIS
2. Picking P-wave arrival times (PhaseNet or STA/LTA)
3. Determining first motion polarity by calculating the cross-correlation values between P-wave and sine wave template

## Project Structure

```
SIN-Polarity/
├── example1-Oklahoma(20141010)/          # Example 1: Oklahoma earthquake
│   ├── 2014-10-10-mb4.5-Oklahoma/        # Raw waveform data (MSEED) and instrument response files (XML)
│   ├── phasenet/                         # PhaseNet source code for P-wave arrival picking
│   ├── model/                            # Pre-trained PhaseNet model
│   ├── SIN_Polarity_Step1_Fetch_data.py               # Step 1: Download data from IRIS
│   ├── SIN_Polarity_Step2_First_arrival_picking.py    # Step 2: Pick P-wave arrival
│   ├── SIN_Polarity_Step3_Polarity_identification.py  # Step 3: Determine polarity
│   └── tool_plot.py                      # Visualization tool
│
├── example2-Myanmar(20250328)/           # Example 2: Myanmar earthquake
│   ├── 2025-03-28-mw7.7-Myanmar/         # Raw waveform data (MSEED) and instrument response files (XML)
│   ├── SIN_Polarity_Step1_Fetch_data.py
│   ├── SIN_Polarity_Step2_First_arrival_picking.py
│   ├── SIN_Polarity_Step3_Polarity_identification.py
│   └── tool_plot.py
│
└── README.md
```

## Example Data

The example folders contain seismic waveform data used in this study:

- **`2014-10-10-mb4.5-Oklahoma/`**: Oklahoma earthquake (2014-10-10, mb4.5)
  - Three-component waveform data in MiniSEED format
  - Station instrument response files in StationXML format
  - Event parameters file (`eventParameters.xml`)

- **`2025-03-28-mw7.7-Myanmar/`**: Myanmar earthquake (2025-03-28, mw7.7)
  - Three-component waveform data in MiniSEED format
  - Station instrument response files in StationXML format
  - Event parameters file (`eventParameters.xml`)

## Key Dependencies

- **`phasenet/`**: Contains the PhaseNet source code (Zhu and Beroza, 2019) for automatic P-wave arrival time picking. This folder is called by `SIN_Polarity_Step2_First_arrival_picking.py` in example1 when using PhaseNet method.

- **`model/`**: Contains the pre-trained PhaseNet model used for arrival time prediction.

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

Run `SIN_Polarity_Step1_Fetch_data.py` to download three-component seismic waveform data and instrument response files for stations within a specified distance range from the epicenter.

**Key parameters to configure**:
- `eventPath`: Directory for storing downloaded data
- `origin_time`: Event origin time (UTC format)
- `lat`, `lon`: Epicenter coordinates (latitude, longitude)
- `minDisInDeg`, `maxDisInDeg`: Station distance range in degrees
- `clientName`: Data service provider (default: "IRIS")

**Output**:
- MiniSEED waveform files (`*.mseed`)
- StationXML instrument response files (`*.xml`)
- Event parameters file (`eventParameters.xml`)

### Step 2: P-wave Arrival Time Picking

Run `SIN_Polarity_Step2_First_arrival_picking.py` to automatically pick P-wave first-arrival times using PhaseNet or traditional STA/LTA algorithm.

**Processing workflow**:
1. Organize three-component data into separate folders (data1, data2, dataZ)
2. Remove linear trend, mean, and bandpass filtering and calculate theoretical P-wave travel time (ak135 model)
3. Extract time window around theoretical arrival
4. Call PhaseNet or STA/LTA for automatic picking

**Output**: `templates.csv` containing P-wave first-arrival times for each station

### Step 3: P-wave First-Motion Polarity Determination

Run `SIN_Polarity_Step3_Polarity_identification.py` to determine P-wave first-motion polarity on the vertical component using sine-wave template cross-correlation.

**Processing workflow**:
1. Read P-wave arrival times and vertical-component waveform data
2. Extract signal and noise windows, calculate SNR
3. Calculate the period of the sine wave template
4. Generate sine-wave template and perform cross-correlation with P-wave
5. Determine polarity

**Input**:
- `templates.csv`: P-wave first-arrival times (from Step 2 or external sources)
- `dataENZ/`: Three-component waveform data in MSEED format

**Output**: `StationinformationPhasenet.csv` / `StationinformationSTALTA.csv` containing P-wave first-motion polarity for each station (1=positive/'U', -1=negative/'D')

### Quick Path: Directly Use Step 3

If you already have P-wave first arrival times from public datasets or other picking methods, you can directly use `SIN_Polarity_Step3_Polarity_identification.py`.

**Prepare your data**:

1. Prepare `templates.csv` with the following columns:
   - `network`: Network code
   - `station`: Station name
   - `onsetP`: P-wave arrival time (seconds from trace start)

   Example:

   ```
   network,station,onsetP
   GS,OK025,41.54
   ```

2. Prepare three-component waveform data in `dataENZ/` folder:
   - Files should be in MSEED format
   - Naming convention: `{network}.{station}.mseed`
   - Each file should contain three traces (three components)

## Requirements

- Python 3.x
- obspy
- numpy
- pandas
- scipy
- tensorflow (for PhaseNet)
- matplotlib

## Contact

Any questions or advices? Please contact:
- huilianma@stu.cdut.edu.cn  (Huilian Ma)
- jianlongyuan@cdut.edu.cn   (Jianlong Yuan)
- j.yu@cdut.edu.cn           (Jiashun Yu)
- daweigao666@163.com        (Dawei Gao)
- 1716136870@qq.com          (Shaojie Zhang)
- 1048087167@qq.com          (Cong Wang)
- 2533511481@qq.com          (Xinran Fan)
