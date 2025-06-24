Motion correction
==================

Why motion correction
-------------------------

.. _unit_localization_label:

Unit localization
--------------------------

The method of localization of units is well described in Boussard, et al. [1]_. Each unit's 3D position is
denoted as {:math:`x`, :math:`y`, :math:`z`}. Channel positions {:math:`x_c`, :math:`y_c`} on the
probe plane adhere to SpikeGLX/Kilosort conventions, with :math:`z_c` set to
0. Using a monopolar current source model, we compute the position of
each unit via

.. math::

    \underset{x,y,z,\alpha}{\operatorname{argmin}} \sum_{c \in
    \mathcal{C}}\left(\operatorname{ptt}_{c}-\frac{\alpha}{\sqrt{\left(x-x_{c}\right)^{2}+\left(y-y_{c}\right)^{2}+z^{2}}}\right)^{2}

.

Here, :math:`\operatorname{ptt}_c` denotes the peak-to-trough amplitude on
channel :math:`c`. :math:`\mathcal{C}` denotes the channel indices of :math:`n` nearest
channels from the peak channel (the channel with maximum peak-to-trough
value). :math:`n` is 20 in this study. :math:`\alpha` represents the overall signal
magnitude in the model. The optimization is done by MATLAB function
lsqcurvefit. The position :math:`\boldsymbol{y}` for each unit is used for
estimating the motion of the probe.

.. [1] Boussard, Julien, Erdem Varol, Hyun Dong Lee, Nishchal Dethe, and Liam Paninski. “Three-Dimensional Spike Localization and Improved Motion Correction for Neuropixels Recordings.” In Advances in Neural Information Processing Systems, 34:22095-105. Curran Associates, Inc., 2021. https://proceedings.neurips.cc/paper/2021/hash/b950ea26ca12daae142bd74dba4427c8-Abstract.html.


.. _motion_estimation_label:

Motion estimation
--------------------------

Probe motion across recording sessions was estimated using matched unit pairs (identified by the clustering algorithm) and their localized spatial positions. Let :math:`N_s` be the total number of sessions, :math:`N_p` the number of matched unit pairs , and :math:`\boldsymbol{p}_s` the probe position for session :math:`s`. For the :math:`i`-th matched pair (:math:`i = 1, \dots, N_p`), let :math:`s_A^{i}` and :math:`s_B^{i}` denote the sessions from which the two matched units originate, with :math:`y_A^i` and :math:`y_B^i` representing their spatial positions along the probe. The probe positions :math:`\boldsymbol{p} = [p_1, \dots, p_{N_s}]` were estimated by solving:

.. math::

    \boldsymbol{p}^* = \underset{\boldsymbol{p}}{\arg\min} \sum_{i=1}^{N_p} [( y_A^i -y_B^i) - (p_{s_A^i} - p_{s_B^i} ) ]^2

.

This minimizes the discrepancy between the relative displacements of matched units and the inferred probe motion across sessions. The optimization was performed using MATLAB's ``fminunc``, and mean-subtracted probe positions (:math:`\boldsymbol{p}^* - \text{mean}(\boldsymbol{p}^*)`) were used for waveform correction to center displacements around a common reference. 

.. _waveform_correction_label:

Waveform correction
--------------------------

We applied Kriging interpolation method for motion correction, adapting
it to interpolate average waveforms instead of raw recordings. Corrected
waveforms :math:`\tilde{\mathbf{W}}` at probe position :math:`v = \{x, y\}` is the
weighted summation of the original waveforms :math:`\mathbf{W}`, weighted by
spatial proximity. The distance matrices between positions
:math:`\boldsymbol{v}_1` and :math:`\boldsymbol{v}_2` were defined as:

.. math::

    \mathbf{D}_x(\boldsymbol{v}_1,\boldsymbol{v}_2) = \lvert
    \boldsymbol{x}_{1}-\boldsymbol{x}^{T}_2 \rvert

and 

.. math::
    \mathbf{D}_y(\boldsymbol{v}_1,\boldsymbol{v}_2) = \lvert
    \boldsymbol{y}_{1}-\boldsymbol{y}^{T}_2 \rvert

,

where :math:`\mathbf{D}_x(\boldsymbol{v}_1,\boldsymbol{v}_2)` and
:math:`\mathbf{D}_y(\boldsymbol{v}_1,\boldsymbol{v}_2)` is the distance matrix
between positions :math:`\boldsymbol{v}_1` and positions :math:`\boldsymbol{v}_2` in
the x and y direction, respectively. Let
:math:`\boldsymbol{v_\mathcal{C}} = \{\boldsymbol{x_\mathcal{C}}, \boldsymbol{y_\mathcal{C}}\}`
denote the positions of the channels :math:`\mathcal{C}`. Then, the corrected
waveform :math:`\tilde{\mathbf{W}}` at probe position :math:`v` is computed via

.. math::
    
    K(\boldsymbol{v}_1,\boldsymbol{v}_2) =
    e^{-\frac{\mathbf{D}_x}{\sigma_x}-\frac{\mathbf{D}_y}{\sigma_y}}

and 

.. math::
    \tilde{\mathbf{W}}(v) =
    K(v,\boldsymbol{v_\mathcal{C}})K(\boldsymbol{v_\mathcal{C}},\boldsymbol{v_\mathcal{C}})^{-1}\mathbf{W}_\mathcal{C}

,

where :math:`K` is a generalized Gaussian kernel. :math:`\sigma_x` and :math:`\sigma_y`
are two parameters controlling the size of the kernel, with
:math:`\sigma_x = 20` and :math:`\sigma_y = 30` in this paper. For cross-session
consistency, all waveforms were aligned to the mean probe position
(reference probe). The corrected waveform :math:`\tilde{\mathbf{W}}` for unit
$i$ on the reference probe is

.. math::

    \tilde{\mathbf{W}^i}(\boldsymbol{v}_\mathcal{C}) = K(\boldsymbol{v_\mathcal{C}} - \{0, \boldsymbol{p}_{
    \boldsymbol{s}_i}
    \},\boldsymbol{v_\mathcal{C}})K(\boldsymbol{v_\mathcal{C}},\boldsymbol{v_\mathcal{C}})^{-1}\mathbf{W}_\mathcal{C}^i

.

Now, we can compute the corrected waveform similarity again as :ref:`before <waveform_similarity_label>`!

Doing motion correction multiple times
----------------------------------------------

1. change of weights of waveform
2. change of probe positions


.. _non_rigid_correction_label:

Non-rigid correction
---------------------------





