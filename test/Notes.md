# Notes on testinng nr-catalog-tools

## SXS waveforms

1. SXS waveforms are tested against [waveformtools](https://gitlab.com/vaishakp/waveformtools)
2. **waveformtools** loads in SpEC data using the package **scri**.
    1. It then stores it in its own `modes_array`. 
    2. Given extrinsic parameters, evaluates the waveform using a native implementaion.

## GT waveforms

1. GT waveforms are tested against LAL NR data injection infrastructure.
2. This uses the SEOBNRv4 ROM which is a part of the [lalsuite-extra]() repository.
3. The ROM file being>300MB, needs to be accessible to the runners used in testing.
    1. It is infeasible to use git lfs since there is a bandwidth limit of 1GB.
    2. Presently, this is being downloaded inside the runner.

## RIT waveforms

1. RIT waveforms are tested against waveformtools.

