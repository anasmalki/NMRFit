from setuptools import setup, find_packages

setup(
    name='nmrfit',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'pathlib',
        'scipy',
        'lmfit',
        'biopython',
        'plotly',
        'pybroom'
        'tqdm'
    ],
)