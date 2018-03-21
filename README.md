# MSssInteg
Tool/steps for using MountainSort as pre-processing of spike snippets recorded by the Plexon MAP system.

NOTE: Working with .plx files requires the Plexon offline SDK which can be found [here](https://plexon.com/wp-content/uploads/2017/08/OmniPlex-and-MAP-Offline-SDK-Bundle_0.zip). Make sure that toolbox is in your path when running this code.

# Overview
1. Use 'ConvertPLXtoMDA.m' to, as the name indicates, convert the data from the Plexon .plx format, into the .mda format used by MountainSort.

2. Take these .mda files and run them through the MountainSort analysis, taking care to follow the instructions outlined on that system.

3. Retrieve the sorted .mda files and then use 'ConvertAnalyzedMDAtoNEX.m' to, as the name indicates, convert these sorted .mda files into .nex files.

4. Import the .nex files into Offline Sorter and then export as .plx.

5. Open the newly created .plx file in Offline Sorter, change  the channel assignment to Tetrode and then assess the output of MountainSort and identify any missed clusters.

