from setuptools import setup, find_packages

setup(
    name='RMSX: Time series RMSF',
    version='0.1.0',
    author='Finn Beruldsen',
    author_email='fpberuld@cougarnet.uh.edu',
    description='A package for RMSX trajectory analysis using MDAnalysis.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/finn2400/rmsx',
    packages=find_packages(),
    include_package_data=True,
    package_data={'rmsx': ['r_scripts/*.R']},
    install_requires=[
        'MDAnalysis',
        'pandas',
        # Include other dependencies, not sure what else i need
    ],
    scripts=['scripts/rmsx_cli.py'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',  # Choose appropriate license
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)

