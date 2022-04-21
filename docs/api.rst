Processing
==========

.. automodule:: S1_NRB.processor
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        main

NRB
----

.. automodule:: S1_NRB.nrb
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        format
        get_datasets
        create_vrt
        create_rgb_vrt
        calc_product_start_stop
        create_data_mask
        create_acq_id_image

ETAD
----

.. automodule:: S1_NRB.etad
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        process

DEM
---

.. automodule:: S1_NRB.dem
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        prepare
        mosaic


Metadata
========

Extraction
----------

.. automodule:: S1_NRB.metadata.extract
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        get_prod_meta
        vec_from_srccoords
        etree_from_sid
        convert_coordinates
        find_in_annotation
        calc_performance_estimates
        extract_pslr_islr
        get_header_size
        meta_dict

XML Parser
----------

.. automodule:: S1_NRB.metadata.xml
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        product_xml
        source_xml
        parse

STAC Parser
-----------

.. automodule:: S1_NRB.metadata.stac
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        product_json
        source_json
        parse

Ancillary Functions
===================

.. automodule:: S1_NRB.ancillary
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        generate_unique_id
        get_max_ext
        set_logging
        log

Tile Extraction
---------------

.. automodule:: S1_NRB.tile_extraction
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        tiles_from_aoi
        extract_tile
        description2dict
        get_tile_dict

Configuration
=============

.. automodule:: S1_NRB.config
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        get_config
        geocode_conf
        gdal_conf
