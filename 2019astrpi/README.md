# 2019astrpi
Characterization of sound-evoked responses in the posterior striatum.

## Generate database
* `database_generate_no_stats.py`: Generate basic database of cells. It creates `astrpi_all_cells_no_stats.h5` (it takes 30s)
* `database_add_cell_locations.py`: Add location of each recorded cell. It creates `astrpi_all_cells_locations.h5` (it takes 20s)
* `database_select_sound_session.py`: Create a database subset that includes only data with sound sessions. It creates `astrpi_cells_with_soundsessions.h5` (it takes 1s)
* `database_remove_bad_cells.py`: Exclude manually selected as low quality (stored in sj_manually_removed.txt). It creates `astrpi_included_cells.h5` (it takes 1s)
* `database_laser_response.py`: Estimate response to laser. It creates `astrpi_laser_response.h5` (it takes 25s)
* `database_tone_responsive.py`: Estimate responses to pure tones. It creates `astrpi_tones_pval.h5` (it takes 7 min)
* `database_am_responsiveness.py`: Estimate responses to AM sounds. It creates `astrpi_am_pval.h5` (it takes 2 min)
* `database_freq_tuning.py`: Estimate frequency tuning and latency (using stimulus timing from taskontrol). It creates `astrpi_freq_tuning.h5` (it takes 2.5 min)
* `database_tone_latency.py`: Estimate response latency (using stimulus timing from sound detector) and onset-to-sustained ratio. It creates `astrpi_freq_tuning_latency.h5` (it takes 1.5 min)
* `database_am_tuning.py`: Estimate tuning to AM rate. It creates `astrpi_am_tuning.h5` (it takes 30s)

## Figure 1
* `figure_photoid.py`

## Figure 2
* `figure_locations_by_responsiveness.py`

## Figure 3
* `figure_tone_responses.py`

## Figure 4
* `figure_response_dynamics.py`

## Figure 5
* `figure_am_responses.py`

