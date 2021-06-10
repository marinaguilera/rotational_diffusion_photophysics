# Experiemtnal conditions of STARSS experriments

## STARSS1
Experiment day: 2020/07/21. 
Sample: Beads.
### Power densities
- 488 power = 6.69uW 5.6mm BA (8mW settings on Cobolt Software)
- 405 power = 58.3uW 5.6mm BA (2mW settings on Cobolt Software)
- 488nm PSF area at FWMH = 3.80E-10 cm2 @ 220nm FWHM
- 405m PSF area at FWMH = 3.30E-10 cm2 @ 205nm FWMH
- 488 power density = 17.6e3W/cm2 (8mW settings on Cobolt Software)
- 405 power density = 176.7e3W/cm2 (2mW settings on Cobolt Software)

### All parameters
                DO0Setting: [10×1 logical]
                DO1Setting: [10×1 logical]
                DO2Setting: [10×1 logical]
                DO4Setting: [10×1 logical]
                DO3Setting: [10×1 logical]
                AO0Setting: [5×1 double]
                AO3Setting: [5×1 double]
                AO4Setting: [5×1 double]
                AO1Setting: [5×1 double]
                AO2Setting: [5×1 double]
        recordTimeFrom_us_: 100
                    points: 8000
             dwellTime_us_: 0.250000000000000
            windowTime_us_: [10×1 double]
                  averages: 5
                      size: [3×1 double]
                offset_um_: [3×1 double]
             pixelSize_um_: [3×1 double]
         calibration_um_V_: [3×1 double]
                 startTime: '07-21-2020_13:54:34,185'
                   endTime: '07-21-2020_14:01:34,687'
      decisionTimeFrom_us_: 0
        decitionTimeTo_us_: 600
    lowerThreshold_counts_: 2
    upperThreshold_counts_: 5000
             smartScanning: 0
             
>> beads(1).data.H.parj.windowTime_us_ =
   1.0e+03 *
   0.50000
   0.00035
   1.10000
   0.50000
   0.00005
                   0
                   0
                   0
                   0
                   0

>> beads(1).data.H.parj.DO0Setting =
   0
   1
   0
   0
   0
   0
   0
   0
   0
   0

>> beads(1).data.H.parj.DO1Setting =
   1
   1
   1
   1
   0
   0
   0
   0
   0
   0


## STARSS2
Experiment day: 2020/10/22.
Sample: beads.
### Power densities
405 nm power: 37.9uW with 5.6 BA, 2mW settings.
488 nm power: 1.861 uW with 5.6 BA, 2mW settings.
488 nm power: 7.44 uW with 5.6 BA, 8mW settings.
405 nm power density: 114.8e3W/cm2, 2mW settings.
488 nm power density: 4.90e3W/cm2, 2mW settings.
488 nm power density: 19.6e3W/cm2, 8mW settings.

### All parameters
>> beadspw2(1).data.H.parj

                DO0Setting: [10×1 logical]
                DO1Setting: [10×1 logical]
                DO2Setting: [10×1 logical]
                DO4Setting: [10×1 logical]
                DO3Setting: [10×1 logical]
                AO0Setting: [5×1 double]
                AO3Setting: [5×1 double]
                AO4Setting: [5×1 double]
                AO1Setting: [5×1 double]
                AO2Setting: [5×1 double]
        recordTimeFrom_us_: 100
                    points: 4000
             dwellTime_us_: 1
            windowTime_us_: [10×1 double]
                  averages: 1
                      size: [3×1 double]
                offset_um_: [3×1 double]
             pixelSize_um_: [3×1 double]
         calibration_um_V_: [3×1 double]
                 startTime: '10-22-2020_19:26:39,539'
                   endTime: '10-22-2020_19:29:24,042'
      decisionTimeFrom_us_: 0
        decitionTimeTo_us_: 600
    lowerThreshold_counts_: 2
    upperThreshold_counts_: 5000
             smartScanning: 0


>> beadspw2(1).data.H.parj.windowTime_us_ =
   1.0e+03 *
   0.50000
   0.00035
   0.50000
   3.10000
   0.00005
   0
   0
   0
   0
   0

>> beadspw2(1).data.H.parj.DO0Setting =
   0
   1
   0
   0
   0
   0
   0
   0
   0
   0

>> beadspw2(1).data.H.parj.DO1Setting =
   1
   1
   0
   1
   0
   0
   0
   0
   0
   0


## STARSS3
Experiment day: 2021/02/16.
Sample: beads.
### Power densities
405 nm power: 2350uW with 5.6 BA, 60mW settings.
488 nm power: 5.32uW with 5.6 BA, 10mW settings.
405 nm power density: 7.12e6W/cm2, 60mW settings.
488 nm power density: 14e3W/cm2, 10mW settings.
### All parameters
                DO0Setting: [10×1 logical]
                DO1Setting: [10×1 logical]
                DO2Setting: [10×1 logical]
                DO4Setting: [10×1 logical]
                DO3Setting: [10×1 logical]
                AO0Setting: [5×1 double]
                AO3Setting: [5×1 double]
                AO4Setting: [5×1 double]
                AO1Setting: [5×1 double]
                AO2Setting: [5×1 double]
        recordTimeFrom_us_: 0
                    points: 800
             dwellTime_us_: 10
            windowTime_us_: [10×1 double]
                  averages: 3
                      size: [3×1 double]
                offset_um_: [3×1 double]
             pixelSize_um_: [3×1 double]
         calibration_um_V_: [3×1 double]
                 startTime: '02-16-2021_18:54:44,588'
                   endTime: '02-16-2021_18:54:54,587'
      decisionTimeFrom_us_: 0
        decitionTimeTo_us_: 600
    lowerThreshold_counts_: 2
    upperThreshold_counts_: 5000
             smartScanning: 0

>> parj.windowTime_us_ = 
   1.0e+03 *
   0.59990
   0.00015
   0.00010
   0.00015
   0.20000
   3.20000
   0.60015
   0.00015
   0.20000
   3.20000

>> parj.DO0Setting =
   0
   1
   0
   1
   0
   0
   0
   1
   0
   0

>> parj.DO1Setting =
   0
   0
   0
   0
   0
   1
   0
   0
   0
   1
