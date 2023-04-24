# Notes on testinng nr-catalog-tools

## SXS waveforms

1. SXS waveforms are tested against [waveformtools](https://gitlab.com/vaishakp/waveformtools)
2. **waveformtools** loads in SpEC data using the package **scri**.
    1. It then stores it in its own `modes_array`. 
    2. Given extrinsic parameters, evaluates the waveform using a native implementaion.

## GT waveforms

1. GT waveforms are tested against LAL NR data injection infrastructure.


## RIT waveforms

1. RIT waveforms are tested against waveformtools.

