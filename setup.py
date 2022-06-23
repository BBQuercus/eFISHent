"""Setup file for pypi package called efishent."""

# TODO write makefile for basic deployment using
# python setup.py sdist
# twine upload dist/latest-version.tar.gz

import textwrap
from setuptools import find_packages
from setuptools import setup

setup(
    # Description
    name="eFISHent",
    version="0.0.1",
    license="MIT",
    description="RNA FISH oligos/probes design tool.",
    long_description_content_type="text/plain",
    long_description=textwrap.dedent("""TDB."""),
    # Installation
    python_requires=">3.9",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        "biopython",
        "gtfparse",
        "luigi",
        "matplotlib",
        "numpy",
        "pandas",
        "pyarrow",
        "pyomo",
        "tqdm",
    ],
    entry_points={
        "console_scripts": [
            "efishent = eFISHent.cli:main",
            "eFISHent = eFISHent.cli:main",
        ]
    },
    # Metadata
    author="Bastian Eichenberger",
    author_email="bastian@eichenbergers.ch",
    url="https://github.com/bbquercus/efishent/",
    project_urls={
        "Documentation": "https://github.com/BBQuercus/efishent/wiki",
        "Changelog": "https://github.com/BBQuercus/efishent/releases",
        "Issue Tracker": "https://github.com/bbquercus/efishent/issues",
    },
    keywords=["biomedical", "bioinformatics", "RNA probe design"],
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Utilities",
    ],
)
