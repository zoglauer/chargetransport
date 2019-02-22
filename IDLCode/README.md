# Original IDL based charge transport simulations

## Requirements

1. Install GDL 0.9.9 or later
2. Create an output directory "idlsave" in this directory

## Usage

1. Go through the Det_ID from 1 to 12 for run_fac_cal.pro file
2. Type 'gdl'
3. Type '@run_fac_cal'

## Looking at the output files with python3 (requires scipy):

```
python3
>>> from scipy.io.idl import readsav
>>> s = readsav("idlsave/NCT14_D10_efield.idlsave",verbose=True)
```
