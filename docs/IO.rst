Input and output
=================

.. contents:: 
    :local:

Input
-------

``spikeInfo.mat``
+++++++++++++++++++

See how to build it :ref:`here <prepare_the_data_label>`.

``settings.json``
+++++++++++++++++++++

See how to prepare it :ref:`here <prepare_the_data_label>`.

See how to adjust the settings :doc:`here <Change_default_settings>`.


Output
-------------

``Output.mat``
++++++++++++++++

See :ref:`here <output_label>` for the details.

.. _motion_output_label:

``Motion.mat``
++++++++++++++++

Contains a struct variable named ``Motion``, which stores the estimated probe motion across sessions. The fields are:

===========================     =============================               =================
Field name                      Type                                        Explanation  
===========================     =============================               =================
``Linear``                      1 x n_session double                        linear term of the motion
``Constant``                    1 x n_session double                        constant term of the motion. If rigid, it is the probe motion.
``LinearScale``                 1 x 1 double                                0.001 by default. It is used to scale the Y position for numerical stability during motion estimation
===========================     =============================               =================

The probe positions across sessions at a certain Y position ``y`` can be computed via

.. code-block:: MATLAB

    probe_positions = Motion.LinearScale*Motion.Linear*y + Motion.Constant;

.

See :ref:`Non-rigid correction <non_rigid_correction_label>` for details.


``Waveforms.mat``
+++++++++++++++++++

Contains a variable named ``waveforms_corrected``, which is a n_unit x n_channel x n_sample 3D array. It stores the corrected waveforms for all units on all the channels. 


``resultIter.mat``
++++++++++++++++++++

Contains a varaible named ``resultIter``, which is a 1 x ``n_iter`` struct array. ``n_iter`` is the times of motion correction. This variable stores the information in each iteration of motion correction. The fields are:

===========================     =============================               =================
Field name                      Type                                        Explanation  
===========================     =============================               =================
``FeatureNames``                1 x n_feature cell array                    features used in this iteration
``Weights``                     1 x n_feature double                        optimized weights in this iteration
``IdxClusters``                 n_unit x 1 double                           clustering result (cluster index for each unit) in this iteration
``Motion``                      1 x 1 struct                                estimated probe position (See :ref:`Motion <motion_output_label>`) in this iteration
===========================     =============================               =================

``SimilarityMatrix.mat``
+++++++++++++++++++++++++++++

Contains two varaibles named ``feature_names_all`` and ``similarity_matrix_all``.

``similarity_matrix_all`` is a n_unit x n_unit x n_feature 3D array, containing all similarity scores for all pairs of units. ``feature_names_all`` is a 1 x n_feature array indexing the features in ``similarity_matrix_all``. Some features are not used in later clustering process, which will be zeros in ``similarity_matrix_all``.

Intermediate files
+++++++++++++++++++++

These files contains ``ClusterIndices.npy``, ``DistanceMatrix.npy``, ``LinkageMatrix.npy`` and ``HDBSCAN_settings.json``. These files are used to communiate the data between MATLAB and Python when running HDBSCAN algorithm. 


