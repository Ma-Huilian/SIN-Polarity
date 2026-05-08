# SIN-Polarity

A tool for determining P-wave first motion polarity based on sine wave template.

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
a_tool_FetchDataFromIRIS.py   # Step 1: Download data from IRIS
        ↓
b_tool_phasenet.py            # Step 2a: Pick P-wave with PhaseNet
  or b_tool_StaLta.py         # Step 2b: Pick P-wave with STA/LTA
        ↓
c_tool_SIN-Polarity.py        # Step 3: Determine polarity
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

## Output

- `templates.csv`: P-wave arrival times for each station
- `StationinformationPhasenet.csv` / `StationinformationSTALTA.csv`: Polarity results (1='U', -1='D')
- `output_FMsln.csv`: Focal mechanism solutions (if applicable)
