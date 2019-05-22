# HiCtool utility code

We provide a series of utility functions included into [HiCtool_utilities.py](/scripts/HiCtool_utilities.py) that allow to work with HiCtool generated data in the Python environment, such as:

- Loading and saving HiCtool contact matrices (either in compressed and tab-separated format).
- Loading DI values or HMM states.
- Loading topological domain coordinates.

These functions will allow you to use HiCtool data for additional analyses, and eventually re-save to file your processed data in order to be plotted using HiCtool. You could even have contact matrices generated with other software, parse them into a numpy matrix format, and then save them using the saving functions in order to be normalized or plotted for example.

First, open your Python or iPython shell and execute the script:
```Python
execfile('HiCtool_utilities.py')
```
Then, call the function you need (see how to use each function in the function documentation inside the script):

- ``save_matrix`` to save a square and symmetric contact matrix (intra-chromosomal or global matrix) in the HiCtool compressed format.
- ``load_matrix`` to load a square and symmetric contact matrix (intra-chromosomal or global matrix) in the HiCtool compressed format.
- ``save_matrix_rectangular`` to save a rectangular contact matrix (inter-chromosomal) in the HiCtool compressed format.
- ``load_matrix_rectangular`` to load a rectangular contact matrix (inter-chromosomal) in the HiCtool compressed format.
- ``save_matrix_tab`` to save any contact matrix in the tab-separated format.
- ``load_matrix_tab`` to load any contact matrix in the tab-separated format.
- ``load_DI_values`` to load DI values.
- ``load_HMM_states`` to load HMM states.
- ``save_topological_domains`` to save topological domains values in tab-separated format.
- ``load_topological_domains`` to load topological domains values in tab-separated format.




