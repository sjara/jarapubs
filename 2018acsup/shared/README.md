# Scripts for reading photoidentification data from Lakunina et al (2020)

In these scripts, you need to modify the variable `dataDir` to point to the location of your data.
* `show_all_waveforms.py`: Load and plot waveforms for all cells.
* `show_example_cells.py`: Show spike rasters for a few cells and the associated spike waveforms.

## Workflow to go from the raw data to these processed data:
1. `database_generation.py`: Create a database with basic info about clusters.
1. `laser_responsiveness.py`: Load spike times and calculate the response to laser for each cell.
1. `save_data_each_cell.py`: Save the spike times and laser onset times for each cell, as well as a simplified dataframe with information for each cell (including spike waveforms).

## Other scripts:
* `review_rasters.py`: Show one raster at a time given the basic celldb.
* `fix_inforecs.py`: Update inforecs to new versions (post 2021)
* `test_data_consistency.py`: Check that npz files match cells in dataframe.

