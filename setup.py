from setuptools import setup, find_packages

setup(
    name='Linkreg',
    version='0.0.4',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'Linkreg=Linkreg.Linkreg:main',  # "script_name=package.module:function"
        ],
    },
    # Add more metadata as needed
)
