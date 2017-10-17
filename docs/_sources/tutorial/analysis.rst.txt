.. _analysis:

Analysis Pipeline
=================

Tutorial of the basic analysis functions.

.. _velocytoloom:

Velocyto Loom
-------------

Let's start with loading the content of the `.loom` file output into an interactive session of python.

.. code-block:: python

   import velocyto as vcy
   vlm = vcy.VelocytoLoom("YourData.loom")

Different steps of analysis can be carried on by simply calling the methods of the `VelocytoLoom`.
New variables, normalization and parameter will be saved and available as attributes of "VelocytoLoom" the object, while the call does not return anything.
For example normalization and log transformation can be performed by calling the `normalize` method:

.. code-block:: python

    vlm.normalize("S", size=True, log=True)
    vlm.S_norm  # contains log normalized

Furthermore "VelocytoLoom" object supports some ready-made plotting functions.
For example, one of the first checks is spliced/unspliced fractions of the dataset can be done by calling:

.. code-block:: python

    vlm.plot_fractions()

The unspliced fraction should be ~10% of all the molecules.

Preliminary Filtering
---------------------
At this point we can perform feature selection and normalization of the data.
In order to obtain better results the preliminary filtering is usually adapted for each dataset.
However, we implemented the method `default_filter_and_norm` as a quick shortcut to get started. 
The method uses some heuristics to set the thresholds to reasonable values considering the size of the dataset.

.. code-block:: python

    vlm.default_filter_and_norm()

Notice that the method supports limited options in comparison to the full API. For a finer tuning of filtering parameters inspect the source code of the method in the `API page <http://velocyto.org/velocyto.py/_modules/velocyto/analysis.html#VelocytoLoom.default_filter_and_norm>`_

Preparation for gamma fit
-------------------------
For the preparation of the gamma fit we smooth the data using a kNN neighbors pooling approach.
kNN neighbors can be calculated directly in gene expression space or reduced PCA space, using either correlation distance or euclidean distance.
The default procedure kNN graph pooling/smoothing is implemented `default_fit_preparation`, finer control can be achieved explicitly calling the `knn_imputation <http://velocyto.org/velocyto.py/fullapi/api_analysis.html#velocyto.analysis.VelocytoLoom.knn_imputation>`_ method.

.. code-block:: python

    vlm.default_fit_preparation()


Gamma fit and extrapolation
---------------------------
To fit gamma to every gene that survived the filtering step we can just call

.. code-block:: python

    vlm.fit_gammas()

The fit can be visualized by calling `plot_phase_portraits` and listing the gene names:

.. code-block:: python

    vlm.plot_phase_portraits(["Igfbpl1", "Pdgfra"])

The extrapolation can be obtained as follows:

.. code-block:: python

    vlm.predict_U()
    vlm.calculate_velocity()
    vlm.calculate_shift(assumption="constant_velocity")
    vlm.extrapolate_cell_at_t(delta_t=1)

Projection of velocity onto embeddings
--------------------------------------
The extrapolated cell state is a vector in expression space (available as the attribute `vlm.Sx_sz_t`).
One of the most convenient way to visualize the extrapolated state is to project it on a low dimensional embedding that appropriately summarizes the variability of the data that is of interest.
The embedding can be calculated with your favorite method or external package as soon as it is saved as an attribute of the `VelocytoLoom` object.
For example, let's use `scikit-learn` TSNE implementation and make it available as `ts` attribute as following:

.. code-block:: python

    from sklearn.manifold import TSNE
    bh_tsne = TSNE()
    vlm.ts = bh_tsne.fit_transform(vlm.pcs[:, :25])

Now we can project on `vlm.ts` by calling `estimate_transition_prob`.

.. warning::
   For big datasets this code can take long time to run! We suggest to run it on multicore machines (since the implementation is fully multithreaded) 

::

    vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts")
    vlm.calculate_embedding_shift(sigma_corr = 0.05)

In case of very big dataset visualizations a good way to summarize the velocity is to visualize it as velocity field calculated on a grid.

::

    vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=300)
    vlm.plot_grid_arrows(scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True}, min_mass=24, angles='xy', scale_units='xy',
                         headaxislength=2.75, headlength=5, headwidth=4.8, quiver_scale=0.47)



