from setuptools import setup, find_packages
from pathlib import Path
import re


def get_version():
    """Get version number from __init__.py"""
    version_file = Path(__file__).resolve().parent / "yleaf" / "__init__.py"
    version_match = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]", version_file.read_text(), re.M
    )
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Failed to find version string")


with open('requirements.txt') as f:
    required = f.read().splitlines()


setup(
    name='Yleaf',
    version=get_version(),
    packages=find_packages(),
    url='https://github.com/genid/Yleaf.git',

    license='MIT',
    author='Bram van Wersch and Diego Montiel Gonzalez',
    author_email='b.vanwersch@erasmusmc.nl',
    description='Tool for human Y-chromosomal phylogenetic analysis and haplogroup inference.',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU License",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=required,
    python_requires=">=3.6",
    entry_points={
        'console_scripts': [
            'Yleaf=yleaf.Yleaf:main',
            'predict_haplogroup=yleaf.predict_haplogroup:main',
        ],
    },
    include_package_data=True,
    # NOTE: shipping data/ recursively. The "data/*" glob used previously
    # only matched top-level entries inside data/, which silently excluded the
    # hg_prediction_tables/, hg19/, hg38/ subtrees and rendered the wheel
    # unusable without a manual data-dir copy. See MANIFEST.in for the sdist
    # equivalent.
    package_data={
        "yleaf": [
            "config.txt",
            "data/hg_prediction_tables/*",
            "data/hg_prediction_tables/major_tables/*",
            "data/hg19/*",
            "data/hg38/*",
        ]
    }
)
