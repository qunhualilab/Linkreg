from setuptools import setup, find_packages

setup(
    name='Linkreg',
    version='0.0.1',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    # Add more metadata as needed
)
