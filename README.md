# ThreePointBendingAnalisys

This program is intended to perform a full standard analysis of a 3-point bending experiment performed in the biomechanics lab. 

- Preparation

> Before using this program, make sure that you have a full set of force-displacement files collected from the Dynamight, 5866, or Electropulse, in a single folder. It helps for organization if each file is named after the bone's experimental ID, such as 106-wt.txt or 106-wt.lvm 

> If you've also performed CT scans to capture the most likely fracture area of the bone, be sure so separate these scans using ContouringGUI so that you have one bone in each folder of DICOM files, named the same way as the text files containing force-displacement data. This way they will be easily matched with the force-displacement data.

- Running the software
> Set the beam width appropriately; 7mm is standard for the Dynamight.

> Select the appropriate testing rig from the dropdown menu - this tells the software which format to look for.

> Load Rig Files. Sorting by filename is the best way if you prepared your files well.

> Load DICOM Files. Select a parent folder containing all the folders which each contain the DICOM files for one bone each.

> A list of rig files and DICOM files will be populated. Make sure the columns match. If needed, supply a map that reorders the second list to the first.

> Click Set Results Path to tell the software where to write results.

> Click Start Analysis to begin analyzing. Select the locations required when prompted. See Michael Brodt or Matt Silva for guidance on how to select these locations if you are unfamiliar.

- Results

> View the results file generated in the set location, or in the Force-displacement folder if writing to the requested folder failed.
