from setuptools import setup, find_packages

setup(
    name='Linkreg',
    version='0.0.3',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'Linkreg=Linkreg.Linkreg:Linkreg',  # "script_name=package.module:function"
        ],
    },
    # Add more metadata as needed
)
