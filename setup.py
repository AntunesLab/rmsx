from setuptools import setup, find_packages

# Read the long description from README.md
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='rmsx-timeseries-rmsf',  # Changed to lowercase and hyphens for compatibility
    version='0.1.0',
    author='Finn Beruldsen',
    author_email='fpberuld@cougarnet.uh.edu',
    description='A package for RMSX trajectory analysis using MDAnalysis.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/finn2400/rmsx',
    packages=find_packages(),
    include_package_data=True,  # Ensures package_data is included
    package_data={
        'rmsx': ['r_scripts/*.R'],  # Includes all .R files in r_scripts
    },
    install_requires=[
        'MDAnalysis>=2.0.0',
        'pandas>=1.1.0',
        'plotly>=4.14.3'  # Moved plotly to install_requires if it's essential
        # Add other essential dependencies here if any
    ],
    scripts=['scripts/rmsx_cli.py'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',  # Ensure this matches your LICENSE file
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)



# from setuptools import setup, find_packages
#
# setup(
#     name='RF',
#     version='0.1.0',
#     author='Finn Beruldsen',
#     author_email='fpberuld@cougarnet.uh.edu',
#     description='A package for RMSX trajectory analysis using MDAnalysis.',
#     long_description=open('README.md').read(),
#     long_description_content_type='text/markdown',
#     url='https://github.com/finn2400/rmsx',
#     packages=find_packages(),
#     include_package_data=True,
#     package_data={'rmsx': ['r_scripts/*.R']},
#     install_requires=[
#         'MDAnalysis',
#         'pandas'
#         # Include other dependencies, i think that should be all of it.,
#     ],
#     extras_require={
#         'viz': [
#             'plotly>=4.14.3'
#         ],
#
#     },
#     scripts=['scripts/rmsx_cli.py'],
#     classifiers=[
#         'Programming Language :: Python :: 3',
#         'License :: OSI Approved :: MIT License',  # Choose appropriate license ...
#         'Operating System :: OS Independent',
#     ],
#     python_requires='>=3.6',
# )
#
