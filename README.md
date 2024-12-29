# Kilomatch

[![View Kilomatch on GitHub](https://img.shields.io/badge/GitHub-Kilomatch-blue.svg)](https://github.com/jiumao2/Kilomatch)

## Installation

- Download the code from the [Kilomatch](https://github.com/jiumao2/Kilomatch)
- The development is based on MATLAB R2022b and the following toolboxes should be installed:
    - `Statistics and Machine Learning Toolbox`
    - `Optimization Toolbox`
    - `Parallel Computing Toolbox`
- Python 3.9 or later with `hdbscan` package installed.

```shell
conda create -n hdbscan python=3.10
conda activate hdbscan
pip install scikit-learn
pip install hdbscan
```

## How to use it

### Prepare the data

- The data should be an 1 x n struct array named `spikeInfo` with the following fields:
    - `SessionIndex`: 1 x 1 int scalar indicating the session. It should start from 1 and be coninuous without any gaps.
    - `SpikeTimes`: 1 x n double array in millisecond.
    - `Waveform`: a n_channel x n_sample matrix of the mean waveform in uV. All units should share the same channels.
    - `Xcoords`: a n_channel x 1 double array of the x coordinates of each channel.
    - `Ycoords`: a n_channel x 1 double array of the y coordinates (depth) of each channel.
    - `Kcoords`: a n_channel x 1 double array of the shank index of each channel (not used so far).
    - `PETH`: recommended but not required, a 1 x n double array of the peri-event time histogram.

- The data should be saved in a `.mat` file and specified in the `settings.json`.
- Edit the `settings.json` file to specify the `path_to_data`, `output_folder` and `path_to_python`.
- Edit the other parameters in the `settings.json` file to suit your data.

### Run the code

- Edit the `path_kilomatch` and `path_settings` in `mainKilomatch.m` and run.

### About the output

- All the temporary files, results, and figures will be saved in the `output_folder` specified in the `settings.json` file.
- The results will be saved in the `output_folder` as `Output.mat`, which will contain the following fields:
    - `NumUnits`: 1 x 1 int scalar of the number of units included in the analysis.
    - `NumSession`: 1 x 1 int scalar of the number of sessions included in the analysis.
    - `Sessions`: 1 x n_unit int array of the session index for each unit.
    - `Params`: a struct of the parameters used in the analysis.
    - `Locations`: a n_unit x 3 double array of the estimated x, y, and z coordinates of each unit.

    - `NumClusters`: 1 x 1 int scalar of the number of clusters found (each cluster has at least 2 units).
    - `IdxCluster`: 1 x n_unit int array of the cluster index for each unit. `IdxCluster = -1` means the unit is not assigned to any cluster.
    - `ClusterMatrix`: a n_unit x n_unit logical matrix of the cluster assignment. `ClusterMatrix(i,j) = 1` means unit `i` and `j` are in the same cluster.
    - `MatchedPairs`: a n_pairs x 2 int matrix of the unit index for each pair of units in the same cluster.  
    - `IdxSort`: a 1 x n_unit int array of the sorted index of the units computed from hierarchical clustering algorithm.

    - `SimilarityNames`: a 1 x n_features cell of the names of the similarity metrics used in the analysis.
    - `SimilarityAll`: a n_pairs x n_features double matrix of the similarity between each pair of units. The pairs can be found in `SimilarityPairs`.
    - `SimilarityPairs`: a n_pairs x 2 int matrix of the unit index for each pair of units.
    - `SimilarityWeights`: a 1 x n_features double array of the weights of the similarity metrics computed from IHDBSCAN algorithm.
    - `SimilarityThreshold`: a 1 x 1 double of the threshold which is used to determine the good matches in `GoodMatchesMatrix` used in auto-curation algorithm.
    - `GoodMatchesMatrix`: a n_unit x n_unit logical matrix of the good matches determined by `SimilarityThreshold`. `GoodMatchesMatrix(i,j) = 1` means unit `i` and `j` is a good match.
    - `SimilarityMatrix`: a n_unit x n_unit double matrix of the weighted sum of the similarity between each pair of units.
    
    - `DistanceMatrix`: a n_unit x n_unit double matrix of the distance between each pair of units.
    - `WaveformSimilarityMatrix`: a n_unit x n_unit double matrix of the waveform similarity between each pair of units.
    - `RawWaveformSimilarityMatrix`: a n_unit x n_unit double matrix of the raw (uncorrected) waveform similarity between each pair of units.
    - `ISI_SimilarityMatrix`: a n_unit x n_unit double matrix of the ISI similarity between each pair of units.
    - `AutoCorrSimilalrityMatrix`: a n_unit x n_unit double matrix of the autocorrelogram similarity between each pair of units.
    - `PETH_SimilarityMatrix`: a n_unit x n_unit double matrix of the PETH similarity between each pair of units.

    - `Motion`: a 1 x n_session double array of the positions of the electrode in each session.
    - `Nblock`: a 1 x 1 int scalar of the number of blocks used in the motion correction.

## Notes

- The project is still under development and fundamental changes may occur.
- Be careful that the waveforms included in this analysis should not be whittened as Kilosort does. Do not use the waveforms extracted from `temp_wh.dat` directly. Do not use `whitening_mat_inv.npy` or `whitening_mat.npy` in Kilosort2.5 / Kilosort3 because they are not what Kilosort used to whitten the data (<https://github.com/cortex-lab/phy/issues/1040>)!
- Please raise an issue if you meet any bugs or have any questions.

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
> DREDge: robust motion correction for high-density extracellular recordings across species. https://www.biorxiv.org/content/10.1101/2023.10.24.563768v1

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

