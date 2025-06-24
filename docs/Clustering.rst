Clustering
============

How to track neurons across days
-----------------------------------

Now we have the methods to compute the similarity between any two units (see :doc:`Features <Features>` and :doc:`Motion correction <Motion_correction>`). Then, how to use it to track neurons across sessions? 

:ref:`motion correction <waveform_correction_label>`

Traditionally, ...

Here, we...

.. _HDBSCAN_label:

HDBSCAN
------------

HDBSCAN is an unsupervised clustering algorithm extending DBSCAN by
converting it into a hierarchical clustering algorithm. We utilized the
Python implementation at https://github.com/scikit-learn-contrib/hdbscan. Parameters in this paper were:
``min_cluster_size`` = 2, ``max_cluster_size`` = maximum session number and
``min_samples`` = 1. The input distance matrix :math:`\mathbf{D}` is defined in
the form:

.. math::

    \begin{aligned}
    \mathbf{D}_{i,j} &={\begin{cases}0&{\text{if }}i=j,\\\frac{1}{1+\tanh(\mathbf{S}_{i,j})}&{\text{else }}\end{cases}}
    \end{aligned}

.

HDBSCAN constructed a single-linkage tree for hierarchical clustering,
which was then optimized using MATLAB function optimalleaforder and used
to sort the units in Fig. 3h-i. UMAP, a dimensional reduction method
preserving the structure, was used to visualize the clustering results.
The precomputed distance matrix :math:`\mathbf{D}` was used as input,
producing a two-dimensional embedding. UMAP parameters were adjusted
(``min_dist`` = 1, ``spread`` = 1, ``n_neighbors`` = 15) for improved visualization.


.. _weight_optimization_label:

Weight optimization
-----------------------

We employed LDA to optimize the weights of individual similarity scores and to define the separation boundary between matched and unmatched unit pairs within the high-dimensional feature space. To limit unmatched pairs, only spatially close unit pairs (within 100 μm in :math:`y` position by default) were included. This analysis was performed using MATLAB's function fitcdiscr. The LDA model assumes that similarity scores for matched and unmatched pairs follow multivariate Gaussian distributions with identical covariance matrices. The model generates a hyperplane that maximizes separation between these two classes. The coefficients of the hyperplane's normal vector served as weights to generate a single, optimized similarity score, reflecting the relative importance of each feature. Projecting data onto this one-dimensional vector maximized the discrimination between matches and non-matches. Additionally, the hyperplane defines a similarity threshold :math:`s_{\text{thres}}`, which is useful in the later curation step.

The initial clustering round identifies matches used for motion correction. To minimize false positives, these matches must also satisfy the LDA decision boundary. The second clustering round generates the final results, followed by the auto curation step.


Iterative clustering algorithm
-----------------------------------











