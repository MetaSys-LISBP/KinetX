# KinetX - High-throughput analysis of (bio)chemical kinetics by NMR

## What is IsoCor?
**KinetX is a scientific software dedicated to the processing of (bio)chemical kinetics experiments by NMR**.
KinetX extracts 1D spectra acquired as pseudo 2D spectra, performs Fourier-transform, phasing and beseline correction on each
1D spectra, and extract chemical shifts, intensity and area of the signal(s) of interest.

It is one of the routine tools that we use at the [MetaSys team](http://www.lisbp.fr/en/research/integrated-metabolism-and-dynamics-of-metabolic-systems.html) and [MetaToul platform](https://www6.toulouse.inra.fr/metatoul_eng/) 
for quantitative analysis of enzymatic or chemical kinetics experiments.

The code is open-source, and available under a GPLv3 license.

## Quick-start
KinetX requires TopSpin 3.0 or higher (tested on 3.0PL3 to 3.5PL7) and run on all plate-forms.

To **install KinetX**:

Unpack the content of KinetX.zip somewhere on your disk, and copy the file 'kinetx.py' in the python program folder of your TopSpin user directory (by default: <TopSpin installation directory>/exp/stan/smr/py/usr)

To **use KinetX**:

KinetX requires as input a pseudo 2D spectra which must be Fourier-transformed in F2.

• Open the pseudo 2D spectra to process.
• Define the window that contains the signal of interest in the 2D spectra (only the most intense peak of the displayed area wil be processed)
• Run KinetX:
· run 'kinetx' in the TopSpin command
· enter the signal name (used for spectra annotation), the total number of experiments to process, and confirm the information provided
• Processing results will be displayed and saved in the res subdirectory of the TopSpin experiments folder


**Results**
The folder res is created in the data subdirectory containing all the experiments, and results are stored in the file 'results.txt':
PeakName	name of the signal
PeakID	ID of the peak in TopSpin
Slice	slice of the spectrum (increment in F2)
F1	F1 chemical shift (ppm)
dwF1	difference of chemical shifts in F1 compared to the reference spectrum, which is the first expno (ppm)
Intensity	peak intensity
resF1	FWHM in F1 dimension (ppm)

Error messages are explicit. After correcting the problem, rerun KinetX.


## Bug and feature requests
If you have an idea on how we could improve KinetX please submit a new *issue*
to [our GitHub issue tracker](https://github.com/MetaSys-LISBP/KinetX/issues).


## Developers guide
### Contributions
Contributions are very welcome! :heart:

Please work on your own fork.

## Reference
Cox, N., et al. An integrated pH meter during reaction monitoring with dual reception 1H, 31P NMR spectroscopy (submitted)

## Authors
Pierre Millard, Neil Cox, Guy Lippens

## Contact
:email: Pierre Millard, millard@insa-toulouse.fr
