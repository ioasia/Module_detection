# Module_detection

Code to run module detection from the proteomics data including: 
1) Network construction
2) Sample specfic module detection
3) Co-regulated module clustering

So far, the code execution order is: 1) ALL2_cor_network_analysis.R, 2) ALL2_run_landscapes.py, 3) ALL2_basins.R.

ALL2_run_landscapes can be run with abundance function on network negated, and this is the default. If you want to run it with non-negated abundance values use:
python ALL2_run_landscapes.py --no-negate

## To do
- Stitch code together to create a workflow
- Add basin visualization part

