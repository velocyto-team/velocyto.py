.. _analysis:

Analysis Pipeline
=================

The following tutorial gives an introduction to the basic analysis functions of the `velocyto` library.

.. _velocytoloom:

Velocyto Loom
-------------

Let's start with loading the content of the `.loom` file  into an interactive session of python.

.. code-block:: python

   import velocyto as vcy
   vlm = vcy.VelocytoLoom("YourData.loom")

Different steps of analysis can be carried on by simply calling the methods of this `VelocytoLoom` object.
New variables, normalized version of the data matrixes and other parameters will be stored as attributes of the "VelocytoLoom" object (method calls will not return any value).
For example normalization and log transformation can be performed by calling the `normalize` method:

.. code-block:: python

    vlm.normalize("S", size=True, log=True)
    vlm.S_norm  # contains log normalized

The :ref:`docstring <analysisapi>` of every function specifies which attributes will be generated or modified at each method call.

"VelocytoLoom" object supports some ready-made plotting functions.
For example, one of the first checks is spliced/unspliced fractions of the dataset can be done by calling:

.. code-block:: python

    vlm.plot_fractions()

You can save the results of your analysis in a serialized object at any time by running:

.. code-block:: python

    vlm.dump_hdf5("my_velocyto_analysis")

In another session you can reload the vlm object by running:

.. code-block:: python

    load_velocyto_hdf5("my_velocyto_analysis.hdf5")

This is similar to what the ``pickle`` module in python standard library is doing but here only the attributes of the ``VelocytoLoom`` object are saved and stored as a hdf5 file.
Notice that the size on disk of the serialized file can change depending on the step of the analysis the object is saved (e.g. pre/post filtering or before/after calculating distance matrixes).

.. note::
    VelocytoLoom object methods operate on the object attributes performing filtering, normalization adn other calcualtion. Therefore the order in which they are run is important to get a meaningful output from ``velocyto``.
    We suggest calling these functions in the order shown in this tutorial or in the :ref:`example notebooks <notebooks>`. 


Start a new analysis - Preliminary Filtering
--------------------------------------------

A good first stem is to clean up the data a bit. 
Let's remove the cells with extremelly low unspliced detection

.. code-block:: python

    vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5))

Let's try now to select relevant features for the downstream analysis.
Let's make velocyto aware of the clusters annotation, if we have some

.. code-block:: python

    vlm.set_clusters(vlm.ca["ClusterName"])

Now using the clustering annotation select the genes that are expressed above a threshold of total number of molecules in any of the clusters.

.. code-block:: python

    vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)
    vlm.filter_genes(by_detection_levels=True)

We can perform feature selection.

.. code-block:: python

    vlm.score_cv_vs_mean(3000, plot=True, max_expr_avg=35)
    vlm.filter_genes(by_cv_vs_mean=True)

Finally we can normalize our data by `size` (total molecule count)

.. code-block:: python

    vlm._normalize_S(relative_size=vlm.S.sum(0),
                 target_size=vlm.S.sum(0).mean())
    vlm._normalize_U(relative_size=vlm.U.sum(0),
                 target_size=vlm.U.sum(0).mean())

For a better understend how to fine tune parameters please consult the `API page <http://velocyto.org/velocyto.py/_modules/velocyto/analysis.html#VelocytoLoom.default_filter_and_norm>`_ or just inspect the docstring of each function.

Preparation for gamma fit
-------------------------
For the preparation of the gamma fit we smooth the data using a kNN neighbors pooling approach.
kNN neighbors can be calculated directly in gene expression space or reduced PCA space, using either correlation distance or euclidean distance.
One example of set of parameters is provided below.

.. code-block:: python

    vlm.perform_PCA()
    vlm.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=3000, b_maxl=1500, n_jobs=16)


Gamma fit and extrapolation
---------------------------
To fit gamma to every gene that survived the filtering step run:

.. code-block:: python

    vlm.fit_gammas()


The fit can be visualized by calling `plot_phase_portraits` and listing the gene names:

.. code-block:: python

    vlm.plot_phase_portraits(["Igfbpl1", "Pdgfra"])

The calcualte velocity and extrapolate the future state of the cells:

.. code-block:: python

    vlm.predict_U()
    vlm.calculate_velocity()
    vlm.calculate_shift(assumption="constant_velocity")
    vlm.extrapolate_cell_at_t(delta_t=1.)

In alternative extrapolation can be performed using the constant unspliced assumption (for more information consult our `preprint <citing>`_)

.. code-block:: python

    vlm.calculate_shift(assumption="constant_unspliced", delta_t=10)
    vlm.extrapolate_cell_at_t(delta_t=1.)

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

    vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,
                                 n_neighbors=3500, knn_random=True, sampled_fraction=0.5)
    vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)

In case of very big dataset visualizations a good way to summarize the velocity is to visualize it as velocity field calculated on a grid.

::

    vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=300)
    plt.figure(None,(20,10))
    vlm.plot_grid_arrows(quiver_scale=0.6,
                        scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True}, min_mass=24, angles='xy', scale_units='xy',
                        headaxislength=2.75, headlength=5, headwidth=4.8, minlength=1.5,
                        plot_random=True, scale_type="absolute")








