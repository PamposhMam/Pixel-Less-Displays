.jpg and .bmp files are only here for testing, you can pop in whatever image you want

Steps to run:
1) Create an edge-only image using 'Edge_Images.bmp'
2) Put the file name into MAY2025_CORRECTIONS_WORKINGONATM.mat'. You will receive a hologram.

Functions:
Arc2: somewhat buggy method of taking 3 coordinates, converting them into a form usable by the SLERP code, and then passing a hologram back into the main script
LineHolo: Line holgoram generator
Simpsons_Rule: Numerical integrator for Arc2
