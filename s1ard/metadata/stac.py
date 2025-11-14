import os
import re
import sys
import shutil
import pystac
from spatialist.ancillary import finder
import logging

log = logging.getLogger('s1ard')


def make_catalog(directory, product_type, recursive=True, silent=False):
    """
    For a given directory of Sentinel-1 ARD products, this function will create a high-level STAC
    :class:`~pystac.catalog.Catalog` object serving as the STAC endpoint and lower-level STAC
    :class:`~pystac.collection.Collection` objects for each subdirectory corresponding to a unique MGRS tile ID.
    
    WARNING: The directory content will be reorganized into subdirectories based on the ARD type and unique MGRS tile
    IDs if this is not yet the case.
    
    Parameters
    ----------
    directory: str
        Path to a directory that contains ARD products.
    product_type: str
        Type of ARD products. Options: 'NRB' or 'ORB'.
    recursive: bool, optional
        Search `directory` recursively? Default is True.
    silent: bool, optional
        Should the output during directory reorganization be suppressed? Default is False.
    
    Returns
    -------
    nrb_catalog: pystac.catalog.Catalog
        STAC Catalog object
    
    Notes
    -----
    The returned STAC Catalog object contains Item asset hrefs that are absolute, whereas the actual on-disk files
    contain relative asset hrefs corresponding to the self-contained Catalog-Type. The returned in-memory STAC Catalog
    object deviates in this regard to ensure compatibility with the stackstac library:
    https://github.com/gjoseph92/stackstac/issues/20
    """
    overwrite = False
    product_type = product_type.upper()
    pattern = fr'^S1[AB]_(IW|EW|S[1-6])_{product_type}__1S(SH|SV|DH|DV|VV|HH|HV|VH)_[0-9]{{8}}T[0-9]{{6}}_[0-9]{{6}}_' \
              fr'[0-9A-F]{{6}}_[0-9A-Z]{{5}}_[0-9A-Z]{{4}}$'
    products = finder(target=directory, matchlist=[pattern], foldermode=2, regex=True, recursive=recursive)
    directory = os.path.join(directory, product_type)
    
    # Check if Catalog already exists
    catalog_path = os.path.join(directory, 'catalog.json')
    if os.path.isfile(catalog_path):
        overwrite = True
        catalog = pystac.Catalog.from_file(catalog_path)
        items = catalog.get_all_items()
        item_ids = [item.id for item in items]
        products_base = [os.path.basename(prod) for prod in products]
        diff = set(products_base) - set(item_ids)
        if len(diff) == 0:
            # See note in docstring - https://github.com/gjoseph92/stackstac/issues/20
            catalog.make_all_asset_hrefs_absolute()
            log.info(f"existing STAC endpoint found: {os.path.join(directory, 'catalog.json')}")
            return catalog
    
    sp_extent = pystac.SpatialExtent([None, None, None, None])
    tmp_extent = pystac.TemporalExtent([None, None])
    
    unique_tiles = list(
        set([re.search(re.compile(r'_[0-9A-Z]{5}_'), prod).group().replace('_', '') for prod in products]))
    products = _reorganize_by_tile(directory=directory, product_type=product_type, products=products,
                                   recursive=recursive, silent=silent)
    
    catalog = pystac.Catalog(id=f'{product_type.lower()}_catalog',
                             description=f'STAC Catalog of Sentinel-1 {product_type} products.',
                             title=f'STAC Catalog of Sentinel-1 {product_type} products.',
                             catalog_type=pystac.CatalogType.SELF_CONTAINED)
    
    for tile in unique_tiles:
        tile_collection = pystac.Collection(id=tile,
                                            description=f'STAC Collection of Sentinel-1 {product_type} products for '
                                                        f'MGRS tile {tile}.',
                                            title=f'STAC Collection of Sentinel-1 {product_type} products for '
                                                  f'MGRS tile {tile}.',
                                            extent=pystac.Extent(sp_extent, tmp_extent),
                                            keywords=['sar', 'backscatter', 'esa', 'copernicus', 'sentinel'],
                                            providers=[pystac.Provider(name='ESA',
                                                                       roles=[pystac.ProviderRole.LICENSOR,
                                                                              pystac.ProviderRole.PRODUCER])])
        catalog.add_child(tile_collection)
        
        items = []
        for prod in products:
            if tile in prod:
                item_path = os.path.join(prod, os.path.basename(prod) + '.json')
                item = pystac.read_file(href=item_path)
                items.append(item)
                tile_collection.add_item(item=item)
            else:
                continue
        
        extent = tile_collection.extent.from_items(items=items)
        tile_collection.extent = extent
    
    # Save Catalog and Collections on disk
    catalog.normalize_and_save(root_href=directory)
    
    # See note in docstring - https://github.com/gjoseph92/stackstac/issues/20
    catalog.make_all_asset_hrefs_absolute()
    
    if overwrite:
        log.info(f"existing STAC endpoint updated: {os.path.join(directory, 'catalog.json')}")
    else:
        log.info(f"new STAC endpoint created: {os.path.join(directory, 'catalog.json')}")
    return catalog


def _reorganize_by_tile(directory, product_type, products=None, recursive=True, silent=False):
    """
    Reorganizes a directory containing Sentinel-1 ARD products based on the ARD type and unique MGRS tile IDs.
    
    Parameters
    ----------
    directory: str
        Path to a directory that contains ARD products.
    product_type: str
        Type of ARD products. Options: 'NRB' or 'ORB'.
    products: list[str] or None, optional
        List of ARD product paths. Will be created from `directory` if not provided.
    recursive: bool, optional
        Search `directory` recursively? Default is True.
    silent: bool, optional
        If False (default), a message for each ARD product is printed if it has been moved to a new location or not.
    
    Returns
    -------
    products_new: list[str]
        An updated list of ARD product paths.
    """
    if products is None:
        parent_dir = os.path.dirname(directory)
        pattern = fr'^S1[AB]_(IW|EW|S[1-6])_{product_type}__1S(SH|SV|DH|DV|VV|HH|HV|VH)_[0-9]{{8}}T[0-9]{{6}}_' \
                  fr'[0-9]{{6}}_[0-9A-F]{{6}}_[0-9A-Z]{{5}}_[0-9A-Z]{{4}}$'
        products = finder(target=parent_dir, matchlist=[pattern], foldermode=2, regex=True, recursive=recursive)
    
    inp = input('WARNING:\n{}\nand the ARD products it contains will be reorganized into subdirectories '
                'based on unique MGRS tile IDs if this directory structure does not yet exist. '
                '\nDo you wish to continue? [yes|no] '.format(directory))
    if inp == 'yes':
        tile_dict = {}
        for prod in products:
            tile = re.search(re.compile(r'_[0-9A-Z]{5}_'), prod).group().replace('_', '')
            if tile in tile_dict and isinstance(tile_dict[tile], list):
                tile_dict[tile].append(prod)
            else:
                tile_dict[tile] = [prod]
        
        tiles = list(tile_dict.keys())
        products_new = []
        for tile in tiles:
            tile_dir = os.path.join(directory, tile)
            os.makedirs(tile_dir, exist_ok=True)
            
            for old_dir in tile_dict[tile]:
                new_dir = os.path.join(tile_dir, os.path.basename(old_dir))
                products_new.append(new_dir)
                
                if os.path.dirname(old_dir) != tile_dir:
                    shutil.move(old_dir, new_dir)
                    if not silent:
                        log.info(f"-> {os.path.basename(old_dir)} moved to {tile_dir}")
                else:
                    if not silent:
                        log.info(f"xx {os.path.basename(old_dir)} already in {tile_dir} (skip!)")
                    continue
        return products_new
    else:
        log.info('abort!')
        sys.exit(0)
