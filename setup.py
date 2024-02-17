from setuptools import setup, find_packages

setup(
    name='Linkreg',
    version='0.1.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'Linkreg=Linkreg.Linkreg:main',  # "script_name=package.module:function"
        ],
    },
    url='https://github.com/qunhualilab/Linkreg',
    # Add more metadata as needed
)
