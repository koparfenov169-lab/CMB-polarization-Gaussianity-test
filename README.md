# Unpolarized-points-CMB


This program creates a list of non-polarized points with their type classification for any given E/B polarization map. The program requires the following inputs:

A file with spin-weighted spherical harmonic coefficients (an accompanying program almfrommap.py is provided to generate such a file from a HEALPix map taken from the Planck Legacy Archive)

A parameter file config.txt with the following parameters (the order of lines and the space after '=' are important):

------------------------------
    alm_filename = alm_E_COM_CMB_IQU-smica_2048_R3.00_full.dat # Name of the file containing spin-weighted spherical harmonic coefficients
    mode = E # Specifies the mode for calculation -- E or B
    original or gaussian? = original # Choice of phases for calculation -- original, or Gaussian while preserving the original power spectrum
    grid size = 2048 # The program uses a uniform grid in theta and phi of size 2 * grid size * ( grid size - 1 ). 
    # It is recommended to set grid size at least 8 times larger than l_max, otherwise grid artifacts may appear
    l_min = 2 # Minimum l for spin-weighted harmonics used in map construction (cannot be less than 2, as the CMB dipole component is undefined)
    l_max = 256 # Maximum l for spin-weighted harmonics used in map construction
    output filename = map.dat # Output filename

The output file contains 3 columns: theta in radians (co-latitude), phi in radians (longitude), point type (1 - saddle, 2 - beak, 3 - comet)
