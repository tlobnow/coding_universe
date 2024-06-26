ELISA README
================
Taylor Lab
02/19/2024

1.  **Libraries and packages**
2.  **Paths to input and output**
    - The **Input** for ELISA analysis is stored in excel files with
      sheets that contain all necessary information for analysis.
    - These excel files are stored in one folder per plate.
    - Each **Plate Folder** is named as Plate_X_YYYYMMDD (X is the plate
      number with default 1, enumerated as needed and the assay date
      YYYYMMDD).
    - Any **Output** (tables/figures) created from the plate information
      is stored in the same folder, including:
      - *Plate_X_YYYYMMDD_Standard_Curve* from the Stimulant Dilution
        series
      - absorption values above ~ 1.3 are excluded from the linear
        regression model, as the machine is optimized for analysis below
        1.3 and would distort the regression trajectory
    - *COMBINED_DATA.csv* = Plate Summary Table
    - *MEANS.csv* = Table with calculated Triplicate Means
    - *MOM_SUBSET.csv* = Table with calculated Mean of Means
    - *STATISTICAL_RESULTS.csv* = Table with calculated p-values and
      significance levels
3.  **Table for Name Conversion**
    - *ELISA_CL_KEY.csv* is a conversion table for cell line names (raw
      data –\> names used on plots)
4.  **Analysis Functions**
    - introduction of all functions written for streamlined analysis
      (filtering, normalization, etc..)
5.  **From raw measurements to estimated sample concentrations**
    - Calculate the mean values per Standard Curve dilution step
    - Plot a Standard curve and fit a linear trend line: This equation
      is then used to estimate IL-2 conc. of our unknown samples
      (**ELISA_Fx()**)
    - Save the plots
6.  **Handling of sample data**
    - Data filtering (if necessary)
    - Normalization (lowest value on plate is subtracted from all
      samples, for relative secretion the values are normalized against
      the highest positive control value)
    - Statistical Analysis (Student’s t-Test is used to assess the
      difference between unstimulated and stimulated condition)
7.  **Plotting of sample data**
