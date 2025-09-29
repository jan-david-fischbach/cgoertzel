from setuptools import setup, Extension

setup(
    setup_requires=["cffi>=1.0.0"],
    packages=['cgoertzel'],
    cffi_modules=["cgoertzel/cgoertzel_build.py:ffibuilder"]
)
