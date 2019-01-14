from setuptools import setup

setup(
    name='gautools',
    packages=['gautools'],
    scripts=[
        'gautools/aml.py',
        'gautools/create_runs.py',
        'gautools/geomRegex.py',
        'gautools/out_to_list.py',
        'gautools/out_to_list_sf.py',
        'gautools/submit_gaussian.py',
        'gautools/xtorun.py',
        'gautools/xyz_to_inp.py',
        'gautools/xyz_to_inpglob.py',
    ],
    url='https://github.com/theavey/QM-calc-scripts',
    license='Apache License 2.0',
    author='Thomas Heavey',
    author_email='thomasjheavey@gmail.com',
    description='A set of scripts that are useful for creating, submitting, '
                'and processing QM calculations',
    install_requires=[
        'MDAnalysis>=0.17.0',
        'thtools',
        'numpy',
        'six',
        'paratemp',
    ],
    classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3',
    ],
    zip_safe=True,
)
