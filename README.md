# Monte-Carlo-Project
Cloud_Plotting.py is the main program to run

    - Accepts arguments D, L, C, N, mat, lmbda, uniform, sample, idx
    - D = cloud fractal dimension
    - L = length of supercube (cm)
    - C = number of cubes each axis is broken into
    - mat = Cloud material, actually only takes "test" for a uniform cloud test, "nopower" or "power" for full non-uniform cloud
    - lmbda = Wavelength of light for test
    - uniform = 0 or 1, uniform density if 1
    - sample = 0 or 1, sample non-uniform density if 1
    - idx = portion of observation points to be tested (used for parallel computing on Longleaf)

For our uniform cloud test, we ran Cloud_Plotting with the arguements: 

2.6 3.086e18 64 32 test 0.55 1 0

For our non-uniform cloud test, we ran Cloud_Plotting with the arguments: 

2.6 3.086e18 64 32 nopower 0.55 0 1 (whichever portion of Aobs 0-79 we were doing)
