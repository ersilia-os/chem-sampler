from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf8") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    install_requires = f.read().splitlines()


setup(
    name="chemsampler",
    version="0.1.0",
    author="Miquel Duran-Frigola",
    author_email="miquel@ersilia.io",
    url="https://github.com/ersilia-os/chem-sampler",
    description="Quickly sample the chemical space around a molecule or a set of molecules",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.8",
    install_requires=install_requires,
    packages=find_packages(exclude=("utilities")),
    entry_points={"console_scripts": []},
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    keywords="generative-models machine-learning chemistry computer-aided-drug-design",
    project_urls={"Source Code": "https://github.com/ersilia-os/chem-sampler"},
    include_package_data=True,
)
