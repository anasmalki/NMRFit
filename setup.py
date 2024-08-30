from setuptools import setup, find_packages

setup(
    name='nmrfit',
    version='0.5',
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
