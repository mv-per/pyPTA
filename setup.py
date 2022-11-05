import os
import sys

from pybind11.setup_helpers import Pybind11Extension
from pybind11.setup_helpers import build_ext
from setuptools import find_packages
from setuptools import setup
from setuptools.command.build_py import build_py as _build_py

cpp_args = ["-std=c++11"]

ext_modules = [
    Pybind11Extension(
        "cPTA",
        sources=["src/pyPTA/cpp_files/module.cpp"],
        include_dirs=["pybind11/include"],
        language="c++",
        extra_compile_args=cpp_args,
    )
]

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


class build_py(_build_py):
    def run(self):
        self.run_command("build_ext")
        return super().run()


setup(
    name="pyPTA",
    version="0.0.1",
    author="Marcus Pereira",
    author_email="mav.pereira@outlook.com",
    description="A python package to calculate PTA mono- and muilticomponent data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    cmdclass={"build_ext": build_ext, "build_py": build_py},
    ext_modules=ext_modules,
    url="https://github.com/mv-per/pyPTA",
    project_urls={
        "Bug Tracker": "https://github.com/mv-per/pyPTA/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.6",
)
