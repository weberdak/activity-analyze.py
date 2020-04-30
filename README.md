# serca_activity
Python script for analyzing SERCA activity assays.

## Example usage

### Bash terminal (Linux)

Copy & paste the following into a text file called analyze.bash:

	# Analyze one assay
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
       --rep3_range 45 60 \
       --rep4_row N \
       --rep4_plate 20191009_plate1.txt \
       --rep4_exclude 0 \
       --rep4_range 40 50 \
       --rep5_row N \
       --rep5_plate 20191009_plate1.txt \
       --rep5_exclude 9 \
       --rep5_range 40 50
       
Run by entering the following into a Bash teminal:
  
  	bash batch.bash
  
Note the replicates 4 and 5 will be ignored ith row set to "N".
  
### Windows Powershell 
  
Same principle as above, but with different line continuation.

Copy & paste the following into a text file (Notepad) called analyze.ps1:

	# Analyze one assay
	python activity_analyze.py `
       --prefix serca_alone `
       --series 8.0 7.5 7.0 6.8 6.6 6.4 6.2 6.0 5.8 5.6 5.4 5.0 `
       --interval 10 `
       --well_serca 0.0044 `
       --well_vol 200 `
       --well_eps 6.22e3 `
       --well_path 0.55 `
       --rep1_row A `
       --rep1_plate 20191009_plate1.txt `
       --rep1_exclude 0 `
       --rep1_range 40 50 `
       --rep2_row B `
       --rep2_plate 20191009_plate1.txt `
       --rep2_exclude 9 `
       --rep2_range 40 50 `
       --rep3_row C `
       --rep3_plate 20191009_plate1.txt `
       --rep3_exclude 0 `
       --rep3_range 45 60 `
       --rep4_row N `
       --rep4_plate 20191009_plate1.txt `
       --rep4_exclude 0 `
       --rep4_range 40 50 `
       --rep5_row N `
       --rep5_plate 20191009_plate1.txt `
       --rep5_exclude 9 `
       --rep5_range 40 50

Open a PowerShell an entering the following:
  
  	.\batch.ps1
  	
  Alternatively, you can right-click on batch.ps1 and "Run with PowerShell"
  
  
## Help menu
  
  	python activity_analyze.py -h
  	
  	Analyze activity assay.

	optional arguments:
  	-h, --help            show this help message and exit
  	--prefix PREFIX       Prefix for output files.
  	--series SERIES [SERIES ...]
                        pCa value from X1 to X12.
  	--interval INTERVAL   Time interval between measurements in seconds.
  	--rep1_row {A,B,C,D,E,F,G,H,N}
                        Row for replicate 1.
  	--rep1_plate REP1_PLATE
                        Plate reader file for replicate 1.
  	--rep1_exclude {0,1,2,3,4,5,6,7,8,9,10,11,12} [{0,1,2,3,4,5,6,7,8,9,10,11,12} ...]
                        Exclude these wells from replicate 1.
 	--rep1_range REP1_RANGE [REP1_RANGE ...]
                        First and last time point to calculate rate for replicate 1.
  	--rep2_row {A,B,C,D,E,F,G,H,N}
                        Row for replicate 2.
  	--rep2_plate REP2_PLATE
                        Plate reader file for replicate 2.
 	--rep2_exclude {0,1,2,3,4,5,6,7,8,9,10,11,12} [{0,1,2,3,4,5,6,7,8,9,10,11,12} ...]
                        Exclude these wells from replicate 2.
  	--rep2_range REP2_RANGE [REP2_RANGE ...]
                        First and last time point to calculate rate for replicate 2.
  	--rep3_row {A,B,C,D,E,F,G,H,N}
                        Row for replicate 3.
  	--rep3_plate REP3_PLATE
                        Plate reader file for replicate 3.
  	--rep3_exclude {0,1,2,3,4,5,6,7,8,9,10,11,12} [{0,1,2,3,4,5,6,7,8,9,10,11,12} ...]
                        Exclude these wells from replicate 3.
  	--rep3_range REP3_RANGE [REP3_RANGE ...]
                        First and last time point to calculate rate for replicate 3.
  	--rep4_row {A,B,C,D,E,F,G,H,N}
                        Row for replicate 4.
  	--rep4_plate REP4_PLATE
                        Plate reader file for replicate 4.
  	--rep4_exclude {0,1,2,3,4,5,6,7,8,9,10,11,12} [{0,1,2,3,4,5,6,7,8,9,10,11,12} ...]
                        Exclude these wells from replicate 4.
  	--rep4_range REP4_RANGE [REP4_RANGE ...]
                        First and last time point to calculate rate for replicate 4.
  	--rep5_row {A,B,C,D,E,F,G,H,N}
                        Row for replicate 5.
  	--rep5_plate REP5_PLATE
                        Plate reader file for replicate 5.
  	--rep5_exclude {0,1,2,3,4,5,6,7,8,9,10,11,12} [{0,1,2,3,4,5,6,7,8,9,10,11,12} ...]
                        Exclude these wells from replicate 5.
  	--rep5_range REP5_RANGE [REP5_RANGE ...]
                        First and last time point to calculate rate for replicate 5.
  	--well_serca WELL_SERCA
                        Concentration of SERCA in well (mg/mL).
  	--well_vol WELL_VOL   Volume of well (uL).
  	--well_eps WELL_EPS   Exinction coefficient of NADH at 340 nm (M-1 cm-1).
  	--well_path WELL_PATH
                        Path length through well (cm).
