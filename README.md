# Kilomatch

[![View Kilomatch on GitHub](https://img.shields.io/badge/GitHub-Kilomatch-blue.svg)](https://github.com/jiumao2/Kilomatch)

## Installation

- Download the code from the [Kilomatch](https://github.com/jiumao2/Kilomatch)
- Python 3.9 or later with `scikit-learn` package installed.

```shell
conda create -n hdbscan python=3.10
conda activate hdbscan
pip install scikit-learn
```

## How to use it

### Prepare the data

- The data should be an 1 x n struct array named `spikeInfo` with the following fields:
    - `SessionIndex`: 1 x 1 int scalar indicating the session. It should start from 1 and be coninuous without any gaps.
    - `SpikeTimes`: 1 x n double array in millisecond.
    - `Waveform`: a n_channel x n_sample matrix of the mean waveform. All units should share the same channels.
    - `Xcoords`: a n_channel x 1 double array of the x coordinates of each channel.
    - `Ycoords`: a n_channel x 1 double array of the y coordinates (depth) of each channel.
    - `Kcoords`: a n_channel x 1 double array of the shank index of each channel (not used so far).

- The data should be saved in a `.mat` file and specified in the `settings.json`.

### Run the code

- Update the `settings.json` file to specify the `path_to_data`, `output_folder` and `path_to_python`.
- Edit the path to `settings.json` in `mainKilomatch.m` and run.

### Check the results

- The results will be saved in the `output_folder` as `Output.mat`.
- The `Output.mat` will contain the following fields:
    - `SessionIndex`: 1 x n double array of the session index.
    - `UnitIndex`: 1 x n double array of the unit index.
    - `ClusterIndex`: 1 x n double array of the cluster index.
    - `Waveform`: a n_channel x n_sample x n_unit matrix of the mean waveform.
    - `Xcoords`: a n_channel x 1 double array of the x coordinates of each channel.
    - `Ycoords`: a n_channel x 1 double array of the y coordinates (depth) of each channel.
    - `Kcoords`: a n_channel x 1 double array of the shank index of each channel (not used so far).

## Notes

- The code is still under development and fundamental changes may occur.
- The parameters may need to be adjusted for different datasets.
- Please raise an issue if you meet any bugs or have any questions.

## References

> [HDBSCAN](https://scikit-learn.org/stable/modules/clustering.html#hdbscan)  
>
> Campello, R.J.G.B., Moulavi, D., Sander, J. (2013). Density-Based Clustering Based on Hierarchical Density Estimates. In: Pei, J., Tseng, V.S., Cao, L., Motoda, H., Xu, G. (eds) Advances in Knowledge Discovery and Data Mining. PAKDD 2013. Lecture Notes in Computer Science(), vol 7819. Springer, Berlin, Heidelberg. Density-Based Clustering Based on Hierarchical Density Estimates  
>
> L. McInnes and J. Healy, (2017). Accelerated Hierarchical Density Based Clustering. In: IEEE International Conference on Data Mining Workshops (ICDMW), 2017, pp. 33-42. Accelerated Hierarchical Density Based Clustering

> [Kilosort](https://github.com/MouseLand/Kilosort)  
> Pachitariu, Marius, Shashwat Sridhar, Jacob Pennington, and Carsen Stringer. “Spike Sorting with Kilosort4.” Nature Methods 21, no. 5 (May 2024): 914–21. https://doi.org/10.1038/s41592-024-02232-7.

> [DREDge](https://github.com/evarol/DREDge)  
> DREDge: robust motion correction for high-density extracellular recordings across species. https://www.biorxiv.org/content/10.1101/2023.10.24.563768v1

> [EasyPlot](https://github.com/jiumao2/EasyPlot)

> [npy-matlab](https://github.com/kwikteam/npy-matlab)

> [JSON+C parsing for MATLAB](https://github.com/seanbone/matlab-json-c/releases/tag/v1.1)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
