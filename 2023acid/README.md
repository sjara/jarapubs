# 2023acid
Effects of DOI on the sound responses of auditory cortical neurons.

# Producing databases
The file `studyparams.py` contains the list of animals used in this study as well as 
relevant file paths and parameters for the database calculations.

* `database_generation.py`: Uses the inforec, ephys data, and behavior data to create a database of cells. Creates `celldb_acid2023_basic.h5`

* `database_cell_locations.py`: Adds columns to the database with location (in CCF coords) of recorded cells. Creates `celldb_acid2023_coords.h5`

# Estimating running
* `estimate_running.py`: estimate running magnitude for each animal. Creates `running_SUBJECT.npz` for each SUBJECT.

# Quantify responses to sounds
* `database_freq_tuning.py`: estimate the responses to pure tones (including frequency tuning). You need to also run this for 'running' and 'notrunning' trials.
** `database_freq_tuning.py running`
** `database_freq_tuning.py notrunning`

# Fig.1
## panel_spike_templates.py: plot templates for example neurons
## database_spikes_snr.py: calculate StdDev and spike amplitudes for included neurons.
## panel_spike_amplitude_histogram.py: plot histogram of spike amplitudes and voltage StdDev.
## panel_spike_traces.py: to plot trace of spikes

# Fig.2
## Cartoon
## figure_overall_firing.py
## python add_cartoons.py 1 pdf

# Fig.3
## figure_freqtuning.py
## python add_cartoons.py 2 pdf

# Fig.4
## figure_oddball_enhancement.py
## python add_cartoons.py 3 pdf

# Fig.5
## figure_modulation_oddball_vs_standard.py

# Fig.6
## figure_early_vs_late_days.py
## extras/generate_early_vs_late_samples.py: Resample 100 times and evaluate OEIsaline>OEIdoi