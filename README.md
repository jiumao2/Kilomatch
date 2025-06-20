# Kilomatch

[![View Kilomatch on GitHub](https://img.shields.io/badge/GitHub-Kilomatch-blue.svg)](https://github.com/jiumao2/Kilomatch)
[![Documentation Status](https://app.readthedocs.org/projects/kilomatch/badge/)](https://kilomatch.readthedocs.io/en/latest/)

A MATLAB toolbox for tracking the neurons across days.

Check out the Python version of Kilomatch [here](https://github.com/jiumao2/pyKilomatch)!  

See the [documentation](https://kilomatch.readthedocs.io/en/latest/) for more details.

## Installation

- Download the code from the [Kilomatch](https://github.com/jiumao2/Kilomatch)
- Development is based on MATLAB R2022b. The following toolboxes are required:
    - `Statistics and Machine Learning Toolbox`
    - `Optimization Toolbox`
    - `Parallel Computing Toolbox`
- Python 3.9 ~ 3.11 with `scikit-learn` and `hdbscan` packages installed.

```shell
conda create -n hdbscan python=3.10
conda activate hdbscan
pip install scikit-learn
pip install hdbscan
```  

Specify the Python executable path in `settings.json` (e.g., `"path_to_python": "path_to_anaconda\\anaconda3\\envs\\hdbscan\\python.exe"`).

## How to use it

### Prepare the data

- The data should be an 1 x n struct array named `spikeInfo` with the following fields:
    - `SessionIndex`: 1 x 1 int scalar indicating the session. It should start from 1 and be coninuous without any gaps.
    - `SpikeTimes`: 1 x n double array in millisecond.
    - `Waveform`: an n_channel x n_sample matrix of the mean waveform in uV. All units should share the same channels.
    - `Xcoords`: an n_channel x 1 double array of the x coordinates of each channel.
    - `Ycoords`: an n_channel x 1 double array of the y coordinates (depth) of each channel.
    - `Kcoords`: an n_channel x 1 double array of the shank index of each channel. It is designed for multi-shank probes such as Neuropixels 2.0.
    - `PETH`: recommended but not required, a 1 x n double array of the peri-event time histogram.

- The data should be saved in a `.mat` file and its path specified in the `settings.json`.
- Edit the `settings.json` file to specify the `path_to_data`, `output_folder` and `path_to_python`.
- Edit the other parameters in the `settings.json` file to suit your data.

### Run the code

- Edit the `path_kilomatch` and `path_settings` in `mainKilomatch.m` or `mainKilomatchMultiShank.m` (for multi-shank probes such as Neuropixels 2.0) and run.

### About the output

- All the temporary files, results, and figures will be saved in the `output_folder` specified in the `settings.json` file.
- For multi-shank probes, each shank will have its own output folder as they were processed individually.  
- The results will be saved in the `output_folder` as `Output.mat`, which will contain the following fields:
    - `NumUnits`: 1 x 1 int scalar of the number of units included in the analysis.
    - `NumSession`: 1 x 1 int scalar of the number of sessions included in the analysis.
    - `Sessions`: 1 x n_unit int array of the session index for each unit.
    - `Params`: a struct of the parameters used in the analysis.
    - `Locations`: an n_unit x 3 double array of the estimated x, y, and z coordinates of each unit.

    - `NumClusters`: 1 x 1 int scalar of the number of clusters found (each cluster has at least 2 units).
    - `IdxCluster`: 1 x n_unit int array of the cluster index for each unit. `IdxCluster = -1` means the unit is not assigned to any cluster.
    - `ClusterMatrix`: an n_unit x n_unit logical matrix of the cluster assignment. `ClusterMatrix(i,j) = 1` means unit `i` and `j` are in the same cluster.
    - `MatchedPairs`: an n_pairs x 2 int matrix of the unit index for each pair of units in the same cluster.  
    - `IdxSort`: a 1 x n_unit int array of the sorted index of the units computed from hierarchical clustering algorithm (single linkage + `optimalleaforder`).

    - `SimilarityNames`: a 1 x n_features cell of the names of the similarity metrics used in the analysis.
    - `SimilarityAll`: an n_pairs x n_features double matrix of the similarity between each pair of units. The pairs can be found in `SimilarityPairs`.
    - `SimilarityPairs`: an n_pairs x 2 int matrix of the unit index for each pair of units.
    - `SimilarityWeights`: a 1 x n_features double array of the weights of the similarity metrics computed from IHDBSCAN algorithm.
    - `SimilarityThreshold`: a 1 x 1 double of the threshold which is used to determine the good matches in `GoodMatchesMatrix` used in the auto-curation algorithm.
    - `GoodMatchesMatrix`: an n_unit x n_unit logical matrix of the good matches determined by `SimilarityThreshold`. `GoodMatchesMatrix(i,j) = 1` means unit `i` and `j` is a good match.
    - `SimilarityMatrix`: an n_unit x n_unit double matrix of the weighted sum of the similarity between each pair of units.

    - `Motion`: a 1 x n_session double array of the positions of the electrode in each session.

    - `RunTime`: a 1 x 1 double of the total run time of the analysis in seconds.

## Notes

- Be careful that the waveforms included in this analysis should not be whitened as Kilosort does. Do not use the waveforms extracted from `temp_wh.dat` directly. Do not use `whitening_mat_inv.npy` or `whitening_mat.npy` in Kilosort2.5 / Kilosort3 because they are not what Kilosort used to whiten the data (<https://github.com/cortex-lab/phy/issues/1040>)!
- Please analyze data from different brain regions like cortex and striatum individually since they might have different drifts and neuronal properties.
- Please raise an issue if you meet any bugs or have any questions. We are looking forward to your feedback!

## References

> [HDBSCAN](https://scikit-learn.org/stable/modules/clustering.html#hdbscan)  
> HDBSCAN - Hierarchical Density-Based Spatial Clustering of Applications with Noise. Performs DBSCAN over varying epsilon values and integrates the result to find a clustering that gives the best stability over epsilon. This allows HDBSCAN to find clusters of varying densities (unlike DBSCAN), and be more robust to parameter selection.
> 
> Campello, R.J.G.B., Moulavi, D., Sander, J. (2013). Density-Based Clustering Based on Hierarchical Density Estimates. In: Pei, J., Tseng, V.S., Cao, L., Motoda, H., Xu, G. (eds) Advances in Knowledge Discovery and Data Mining. PAKDD 2013. Lecture Notes in Computer Science(), vol 7819. Springer, Berlin, Heidelberg. Density-Based Clustering Based on Hierarchical Density Estimates  
>
> L. McInnes and J. Healy, (2017). Accelerated Hierarchical Density Based Clustering. In: IEEE International Conference on Data Mining Workshops (ICDMW), 2017, pp. 33-42. Accelerated Hierarchical Density Based Clustering

> [Kilosort](https://github.com/MouseLand/Kilosort)  
> Fast spike sorting with drift correction  
> 
> Pachitariu, Marius, Shashwat Sridhar, Jacob Pennington, and Carsen Stringer. “Spike Sorting with Kilosort4.” Nature Methods 21, no. 5 (May 2024): 914–21. https://doi.org/10.1038/s41592-024-02232-7.

> [DREDge](https://github.com/evarol/DREDge)  
> Robust online multiband drift estimation in electrophysiology data  
> 
> Windolf, Charlie, Han Yu, Angelique C. Paulk, Domokos Meszéna, William Muñoz, Julien Boussard, Richard Hardstone, et al. “DREDge: Robust Motion Correction for High-Density Extracellular Recordings across Species.” Nature Methods, March 6, 2025. https://doi.org/10.1038/s41592-025-02614-5.

> [EasyPlot](https://github.com/jiumao2/EasyPlot)  
> A MATLAB package for making scientific figures easily

> [npy-matlab](https://github.com/kwikteam/npy-matlab)  
> Experimental code to read/write NumPy .NPY files in MATLAB

> [JSON+C parsing for MATLAB](https://github.com/seanbone/matlab-json-c/releases/tag/v1.1)  
> A simple parser for JSON with Comments written in MATLAB

> [MatlabProgressBar](https://github.com/JAAdrian/MatlabProgressBar)  
> This MATLAB class provides a smart progress bar like tqdm in the command window and is optimized for progress information in simple iterations or large frameworks with full support of parallel parfor loops provided by the MATLAB Parallel Computing Toolbox.  

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

