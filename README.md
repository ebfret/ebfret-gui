ebfret-gui
==

This repository contains a reorganized version of the ebfret-legacy code, that has been designed to interact with a graphical user interface (GUI). A manual with installation and usage instructions is located at

    documentation/ebfret_user_guide.pdf 

This project is open source software. Please report any bugs or feature requests via 

    https://github.com/ebfret/ebfret-gui/issues 

Quick Start Guide
==

Installing ebFRET
--

1.  (Windows Only) Download and install the Microsoft Visual C++ 2010 Redistributable Package (x64) from Microsoft (this is required to get the ebFRET code to run).

    http://www.microsoft.com/en-us/download/confirmation.aspx?id=14632

2.  Download ebFRET from

    https://github.com/ebfret/ebfret-gui/archive/1.0.zip    

3.  Unzipping will create a folder 'ebfret-gui-master'.

4.  Start up Matlab and add ebfret to your Matlab path by going to 'File -> Set Path…', or by typing:

    addpath(genpath('/PATH/TO/ebfret-gui-master/src'))

    Substitute /PATH/TO/ with the folder where the zip file was extracted. Note: sometimes there will be another folder called 'ebfret-gui-master' inside the 'ebfret-gui-master' folder. In this case the directory above should be '/PATH/TO/ebfret-gui-master/ebfret-gui-master/src'

5.  Start the GUI with 

    ebf = ebFRET()

6.  Optionally, to use multiple processors, you can type the following

    matlabpool('open', num_cores)

    where num_cores is the number of processors you would like to use (typically 4)


Running ebFRET
--

1.  Load some datafiles by pressing 'File -> Load' and selection the 'Raw Donor Acceptor Time Series (*.dat)' option. You can load multiple files at once. These files are assumed to be part of one 'group' (e.g. one set of experiments performed under the same conditions). 

2.  Remove photo-bleaching. There is an auto-photobleaching removal algorithm which detects the last drop in either donor or acceptor. There is also a manual setting that can filter by setting a threshold on the donor, acceptor, combined intensity or FRET efficiency.

3.  Now look at the histogram in the bottom left. If there are no outliers then you should see a smooth set of peaks over a range of 0 to 1. If you still have outliers, you have two options:

    a.  Click through the traces using the 'Next' and 'Prev' buttons under the 'Time Series' panel. You can crop each trace by adjusting the 'Min' and 'Max' values. If a trace shows a photo-blinking event or should probably be ignored for some other reason, you can click 'Exclude'

    b.  You can clip outlier points, and automatically exclude traces that contain too many outlier points by selecting 'Analysis -> Clip Outliers' 

4. One you're ready to run, set the min and max number of states in the box on the bottom left. And start analysis by pressing the 'Run' button.

    Analysis can be stopped at any time. When restarted it begins with the lowest number of states. If you'd like to just run the analysis for a specific number of states then select 'States: Current' in the analysis panel on the bottom right. 

5. Output can be saved and re-loaded the ebFRET gui, there are also functions that export viterbi paths and summary statistics to an excel-readable CSV format.

