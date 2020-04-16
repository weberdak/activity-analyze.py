# serca_activity
Python script for analyzing SERCA activity assays.

## Example usage

	python activity_analyze.py \
       --prefix serca_alone \
       --series 8.0 7.5 7.0 6.8 6.6 6.4 6.2 6.0 5.8 5.6 5.4 5.0 \
       --interval 10 \
       --well_serca 0.0044 \
       --well_vol 200 \
       --well_eps 6.22e3 \
       --well_path 0.55 \
       --rep1_row A \
       --rep1_plate 20191009_plate1.txt \
       --rep1_exclude 0 \
       --rep1_range 40 50 \
       --rep2_row B \
       --rep2_plate 20191009_plate1.txt \
       --rep2_exclude 9 \
       --rep2_range 40 50 \
       --rep3_row C \
       --rep3_plate 20191009_plate1.txt \
       --rep3_exclude 0 \
       --rep3_range 45 60