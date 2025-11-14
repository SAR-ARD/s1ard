Processing
----------

.. automodule:: s1ard.processor
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        main

SAR
^^^
`s1ard` offers a mechanism to plug in different SAR processors.
The software offers a `snap` reference implementation module that can be translated to other processors.
All "main interface" functions need to be implemented so that `s1ard` can fully interact with the module.

SNAP
++++

.. automodule:: s1ard.snap
    :members:
    :undoc-members:
    :show-inheritance:

    .. rubric:: main interface

    .. autosummary::
        :nosignatures:

        config_to_string
        get_config_keys
        get_config_section
        get_metadata
        process
        translate_annotation


    .. rubric:: processor-specific functions

    .. autosummary::
        :nosignatures:

        grd_buffer
        look_direction
        nrt_slice_num
        pre

ARD
^^^

.. automodule:: s1ard.ard
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        calc_product_start_stop
        create_acq_id_image
        create_data_mask
        create_rgb_vrt
        create_vrt
        format
        get_datasets
        wind_normalization

ETAD
^^^^

.. automodule:: s1ard.etad
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        process

OCN
^^^

.. automodule:: s1ard.ocn
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        extract
        gapfill