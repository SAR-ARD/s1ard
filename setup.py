from setuptools import setup, find_namespace_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='S1_NRB',
    setup_requires=['setuptools_scm'],
    use_scm_version=True,
    description="A prototype processor for Sentinel-1 Analysis Ready Data (ARD) backscatter products",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/SAR-ARD/S1_NRB",
    author="John Truckenbrodt, Marco Wolsza",
    author_email="john.truckenbrodt@dlr.de",
    packages=find_namespace_packages(where='.'),
    include_package_data=True,
    install_requires=['python-dateutil',
                      'gdal>=3.5.0',
                      'click>=7.1.0',
                      'lxml',
                      'pystac>=1.8.0',
                      'pyroSAR>=0.25.0',
                      'spatialist>=0.12.0',
                      'scipy',
                      's1etad>=0.5.3',
                      's1etad_tools>=0.8.1',
                      'pystac-client>=0.7.0',
                      'numpy',
                      'asf_search',
                      'pyproj'],
    extras_require={
          'docs': ['sphinx', 'sphinxcontrib-bibtex', 'nbsphinx', 'sphinx_rtd_theme', 'sphinx-toolbox', 'ipython'],
          'tests': ['pytest', 'requests']
    },
    python_requires='>=3.8',
    license='MIT',
    zip_safe=False,
    entry_points={
        'console_scripts': ['s1_nrb=S1_NRB.cli:cli']
    }
)
