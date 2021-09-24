from setuptools import setup, find_packages
from dvoc_model import __version__ as version

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="dvoc_model",
    version=version,
    author="Donny Zimmanck",
    author_email="dzimmanck@enphaseenergy.com",
    description="Simple Python library for modeling dispatch-able virtual oscillators",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://https://github.com/dzimmanck/dvoc_model",
    packages=find_packages(),
    install_requires=[],
    python_requires=">=3.6",
)