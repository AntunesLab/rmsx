from pathlib import Path
from setuptools import setup, find_packages

README = Path(__file__).with_name("README.md").read_text(encoding="utf-8")

setup(
    name="rmsx-and-flipbook",               # final package name on PyPI/pip
    version="0.1.1",                        # bump to flush old metadata/caches
    description="RMSX trajectory analysis + FlipBook (MDAnalysis + R + ChimeraX)",
    long_description=README,
    long_description_content_type="text/markdown",
    author="Finn Beruldsen",
    author_email="fpberuld@cougarnet.uh.edu",
    url="https://github.com/AntunesLab/rmsx",
    license="MIT",
    packages=find_packages(exclude=("tests*", "docs*", "examples*")),
    include_package_data=True,
    package_data={"rmsx": ["r_scripts/*.R"]},  # ship the R helpers
    install_requires=[
        "MDAnalysis>=2.0.0",
        "pandas>=1.1.0",
        "plotly>=4.14.3",
    ],
    scripts=["scripts/rmsx_cli.py"],        # simplest working CLI hookup
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
