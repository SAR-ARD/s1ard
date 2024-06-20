API Documentation
=================

Configuration
-------------

.. automodule:: s1ard.config
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        gdal_conf
        snap_conf
        get_config

Processing
----------

.. automodule:: s1ard.processor
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        main

SNAP
^^^^

.. automodule:: s1ard.snap
    :members:
    :undoc-members:
    :show-inheritance:

    .. rubric:: core processing

    .. autosummary::
        :nosignatures:

        process
        geo
        grd_buffer
        gsr
        mli
        pre
        rtc
        sgr
        look_direction

    .. rubric:: ancillary functions

    .. autosummary::
        :nosignatures:

        find_datasets
        get_metadata
        postprocess
        nrt_slice_num

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

DEM
^^^

.. automodule:: s1ard.dem
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        mosaic
        prepare

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

Tile Extraction
---------------

.. automodule:: s1ard.tile_extraction
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        aoi_from_scene
        aoi_from_tile
        description2dict
        tile_from_aoi

Ancillary Functions
-------------------

.. automodule:: s1ard.ancillary
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        buffer_min_overlap
        buffer_time
        check_scene_consistency
        check_spacing
        compute_hash
        generate_unique_id
        get_max_ext
        group_by_time
        set_logging
        vrt_add_overviews

Scene Search
------------

.. automodule:: s1ard.search
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        ASF
        ASFArchive
        STACArchive
        asf_select
        check_acquisition_completeness
        collect_neighbors
        scene_select

Metadata
--------

Extraction
^^^^^^^^^^

.. automodule:: s1ard.metadata.extract
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        calc_enl
        calc_geolocation_accuracy
        calc_performance_estimates
        calc_pslr_islr
        copy_src_meta
        find_in_annotation
        geometry_from_vec
        get_header_size
        get_prod_meta
        get_src_meta
        meta_dict

XML
^^^

.. automodule:: s1ard.metadata.xml
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        parse
        product_xml
        source_xml

STAC
^^^^

.. automodule:: s1ard.metadata.stac
    :members:
    :undoc-members:
    :show-inheritance:

    .. autosummary::
        :nosignatures:

        parse
        product_json
        source_json
        make_catalog
