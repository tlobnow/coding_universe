ELISA_BUNDLE
================
Finn Lobnow
, Last edited on 27 November 2024

- [Prerequisites](#prerequisites)
- [Example Files](#example-files)
- [What you can find in the ELISA_BUNDLE
  folder](#what-you-can-find-in-the-elisa_bundle-folder)
- [Analysis Steps](#analysis-steps)
  - [How to prepare your excel sheets for
    Analysis](#how-to-prepare-your-excel-sheets-for-analysis)
  - [How to run the example data](#how-to-run-the-example-data)
  - [How to prepare and run your own
    data](#how-to-prepare-and-run-your-own-data)
  - [Plot and save your data](#plot-and-save-your-data)
  - [Troubleshooting](#troubleshooting)

# Prerequisites

This part should be interesting to you once you have performed an ELISA
as described in . You should have data that ideally contains a positive
and negative control, as well as a plate blank, and you should have the
raw measurements from the plate reader.

# Example Files

If you have access to data-tay, you can find a template file for the
analysis right here: /Volumes/TAYLOR-LAB/ELISA_Analysis_Resources/ This
folder contains the ELISA_BUNDLE with a complete example set. Transfer
this folder to your computer (drag and drop or copy and paste to a
location of your choice).

# What you can find in the ELISA_BUNDLE folder

- archive: contains a static example output folder
- input: put your prepared excel files here
- example_data: I have prepared an example data set to ensure that you
  can run the script without any errors
- ELISA_CL_KEY.csv
  - a static copy of the cell line information
  - the most up-to-date version is stored on my GitHub
  - if you need to make changes, you can add your new cell lines to this
    file
  - column documentation:
    - CELL_LINE: what is the name you provided on the raw plate file?
    - CL_NUMBER: what is the number in our cell line database?
    - CL_NAME_ON_PLOT: what name should be displayed on the plot?
    - PLOT_ID: concatenated CL_NAME_ON_PLOT and EXPRESSION_LVL
    - PURPOSE: what was the purpose of the cell line? Was it a sample or
      a control?
    - PROMOTER: was a strong or weak promoter used? e.g. dSV40 or SFFV
    - ORIGIN: are the components from Mouse or something else?
    - PLOTTING_COLOR: what color should the cell line have in the plot?
    - EXPRESSION_LVL: were the cells sorted for endogenous expression or
      high expression?
    - ORDER_NO: how do you want the cell lines to be ordered in the
      plot?
    - INITIALS: who generated the cell line?
    - INFO: anything noteworthy?
    - COHORT: handy if you also analyzed this cell line with our [Image
      Analysis
      Pipeline](https://github.com/MJ-Taylor-Lab/dynamics-pipeline)
- ELISA_BUNDLE.Rmd
  - a copy of what you are reading here
- ELISA_TEMPLATE.R
  - the script that will be run to analyze your data
  - keep this template and copy the folder for each analysis to keep a
    record of your analysis

# Analysis Steps

## How to prepare your excel sheets for Analysis

Take a look at the example data in the **example_data** folder. You can
use this as a template for your own data. Please make sure that you name
the files and sheets as described in the example data.
(YYYYMMDD_Plate_X) The script will look for the following sheets in your
excel files:

- MEASUREMENTS: contains the raw measured values per plate
  - each plate should contain at least two standard replicates, ideally
    also a positive + negative control
- CELL_LINES: contains the cell line names, standard curve
  concentrations, and where the BLANK is located
- STIM_DAYS: contains the stimulation days (numeric, e.g. 1,2,3,..)
- STIM_TIME: how long you stimulated the cells (e.g. 24 hours)
- STIM_CONCENTRATION: how high was the final stimulation concentration
  with your stimulant (IL-1ß, LPS)?
  - normal IL-1ß concentration is 5ng/mL
    - How is that calculated? We dilute 2µL IL-1ß in 10mL RPMI medium
      (essentially 20ng/mL), but since we add 50µL of this IL-1ß
      containing medium to 150µL medium with settled cells, the final
      concentration is 5ng/mL. For LPS stimulations, 100ng/mL were used.
    - Only the stimulated wells are marked, the unstimulated and the
      standards can be left blank or as NA
- CONDITIONS: which wells are stimulated, unstimulated, or used for
  calibration?
  - Stimulated wells –\> STIM
  - Unstimulated wells –\> UNSTIM
  - Wells that were used for calibration –\> CALIBRATION
  - empty wells can be left blank or as NA
- DILUTIONS: did you dilute the samples? (e.g. 1:5 would mean that each
  cell contains a 5)
  - The standard protocol suggests a 1:5 dilution for all sample wells
    (control and other cell lines), but adjust for other dilutions as
    needed
  - Note for ELISA assay planning:
    - dilute the positive control 1:10 to allow all other samples more
      developing time
    - you need to stop the assay when the color reaches a high intensity
    - Measured values ≤ 1.3 are good, anything above should be taken
      with a grain of salt
- PATHWAY: what pathway are you interested in? (e.g. IL-1R, TLR4, etc.)
- STIMULANT: what stimulant did you use? (e.g. IL-1, LPS, etc.)
- Info: any additional information you want to add to the plot

## How to run the example data

You don’t need to change anything in the ELISA_TEMPLATE.R to run the
example data. I woud advise to open the script in R and take a look at
the script to understand what is happening. To run it, just open
ELISA_TEMPLATE.R and click **Source** on the top right of the script
window. This will run the script and create the output folder with the
plots and tables. The output folder should now contain:

- ELISA_rltv.svg: a plot of the relative secretion values that was
  normalized based on the negative control (lowest value of the negative
  control was subtracted) and the positive control (~1)
- ELISA_real.svg: a plot of the real secretion values that were
  extrapolated based on the standard curve
- ELISA_fold_change.svg: a plot of the fold change values that were
  calculated based on the fold change between the stimulated and
  unstimulated wells per cell line. The statistice were performed based
  on the fold change value of the respective negative control.
- ELISA_ANALYSIS.xlsx: a table with the raw data, the normalized values,
  the real values, and the fold change values

## How to prepare and run your own data

To run your own data, you should copy the ELISA_BUNDLE folder and rename
it to your liking. Open the ELISA_TEMPLATE.R script and adjust file
paths to your own data:

- line 60: instead of the example_data folder, you should point to the
  **input** folder that contains your data

Adjust the actual data:

- line 85 (part 1 ADD YOUR DATA & NORMALIZE IT): instead of the example
  data, you prepare your own data (e.g. you could filter for the plate
  dates of interest, stimulation days, whatever is necessary for you).
  You should also provide the respective positive and negative control
  names that you used in your experiment. This is important for the
  normalization of the data. Make sure you normalize the different data
  sets you provided.
- line 115 (part 2 BIND THE DATA OF INTEREST): bind the normalized data
  batches you wish to visualize below

## Plot and save your data

- If you run the **RELATIVE_SECRETION** or the **REAL_SECRETION**
  chunks, you can plot the data as you wish
  - I have provided the simplest use case with a easy-plot function, but
    there is also a ggplot chunk that you can use to plot the data to
    your needs
- If you want to change the colors of the cell lines, you can do so in
  the **ELISA_CL_KEY.csv** file
- If you run the **SAVE** chunk, you can see that the plots and a new
  Excel file are created in the output folder
- If you run the **FOLD_CHANGE** chunk, you can see that the fold change
  plot is created in the output folder as well

## Troubleshooting

If you run into problems, please make sure that you have the correct
file paths and that the data is correctly formatted. If you still run
into problems, please reach out to me (<flobnow@gmail.com>) and I will
help you out.
