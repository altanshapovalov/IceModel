**READ FIRST:** here you can find brief descriptions of each script in the repository. Instead of uploading the entire working repo that I used during my masters with IceModel, I chose a dozen scripts that best represent my most up to date coding practices and abilities.

**1. a_McdAlbExtract.m** – read .csv files extracted from Google Earth Engine, which contain albedo data from MCD43A3 (NASA’s MODIS). 1) Read csv file, 2) identify any gaps in the data [if gaps are small – patch], and 3) organize data files so that they can be inserted into the IceModel.

**2. compare_mcd43a3_singlepixel.m** – compare albedo from datasets with different spatial resolution. 1) Import data for 5 km and 10 km, 2) select a desired time period [e.g., summer], 3) compare using scatter plots, r-squared, etc., and 4) produce comprehensive plots.

**3. fig_x3_partition.m** – produce final figures for thesis with most in-depth analysis compared to other figures. 1) Establish setting [e.g., select which station’s data will be analyzed, note albedo thresholds, etc.], 2) import melt and albedo data, 3) identify when albedo data is above/below set thresholds, 4) organize that divided data into five distinct periods, 5) display the partitioned information in a comprehensive table for albedo, melt, and time [i.e., duration], 6) calculate and display RMSE information, 7) produce RMSE bar graphs with flexible switches for albedo/melt, and 8) generate box plots with flexible switches for albedo/melt/duration. P.S.: best suited for analyzing+plotting more than one automated weather station data.

**4. fig_x3_upd_oct2024.m** – same as “3. fig_x3_partition.m” but tailored to show more detailed information for a single automated weather station. More specifically, the outermost loop in “3. fig_x3_partition.m” focuses on the weather stations, while in this script it revolves around the available years of data for a given station.

**5. g_x_comparing_tair.m** – compare air temperature (2 m) data of different models relative to data collected by automated weather stations. 1) Import temperature data, 2) extract data for a specific time frame, 3) have the option to produce yearly time series for visual inspection, 4) generate average time series and plot together all datasets, 5) calculate RMSE, and 6) make scatter plots and display as a 3x3 panel.

**6. new_aws_check.m** – read, inspect, and save meteorological data from a dozen automated weather stations on the Greenland Ice Sheet (PROMICE/GC-Net). 1) Read nc-files obtained from PROMICE, 2) select correct attributes for further analysis, 3) apply a series of organizational measures to keep the data files ready [e.g., omit particular columns, format time, convert units, update variable properties, etc.], 4) separate data into individual years, 5) perform an in depth analysis of number of data gaps in each year, 6) patch gaps that are small enough, 7) flag and discard years with too many gaps by displaying it in the code and with generated figures, and 8) plot each variable for every year and weather station for visual inspection.

**7. preprocess_melt_newaws.m** – reorganize meltwater data produced by the IceModel and prepare it for further analysis and plotting. 1) Import IceModel’s melt data, 2) ensure if data transfer was successful by performing value comparison between files, 3) save average, 25th, and 75th percentile data, 4) have the option to plot and inspect the percentile data, and 5) export data with unique and comprehensive names.

**8. read_MAR.m** – read, organize, and save albedo data from the regional climate model MAR (Xavier Fettweis). 1) Read .nc files, 2) select correct coordinates from a 4D variable, 3) expand daily data into hourly, and 4) organize data files for each year so that it is ready to be inserted into IceModel.

**9. read_MERRA2_glc.m** – similar to “8. read_MAR.m” but adjusted to read .nc files of MERRA-2’s M2T3NXGLC product (NASA) to extract albedo, surface temperature (2m), and shortwave radiation.

**10. read_Promice_AWS.m** – similar to “8. read_MAR.m” in the sense that this script also reads .nc files and prepares albedo data for insertion into IceModel; however, this shows a more detailed approach to analyzing multiple weather stations and performs more work in identifying and patching data gaps, which are more common in physical instruments (i.e., weather stations instruments such as pyranometers) rather than model outputs (e.g., data from MAR or MERRA-2).

**11. rename_input_files.m** – an example of a small script that was written for reorganizational and optimization purposes. 1) Import albedo data, 2) rename, 3) ensure that values were unaffected, and 4) export in a repository that IceModel could utilize during simulations.

**12. renaming_ncfiles.m** – rename thousands of .nc files on computer to ensure proper order for reading them and ensure that no data was affected.
