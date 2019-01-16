CLI Logic
=========

The following page describes the API of the logic module.

.. warning::
    The logic module has been introduced as an experimetal feature in version 0.12.0, we forecast that its API and design will chnge considerably in future varsions. 

The logic module offers the possibility to easily change the logic used to count "spliced", "unspliced", "ambiguous" (and other categories) of molecules, without requiring a complete rewrite of the all pipeline.
This might be desired to adapt velocyto to different technologies or to introduce new rules and heuristics to make the counting more accurate or sensitive.

velocyto\.logic module
----------------------

.. _logicapi:

.. automodule:: velocyto.logic
    :members:
    :undoc-members:
    :show-inheritance:
