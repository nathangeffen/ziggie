import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ziggie",
    version="0.1.0",
    author="Nathan Geffen",
    author_email="nathan@simhub.online",
    description="Infectious disease modelling package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nathangeffen/ziggie",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",

        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",

        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",

        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3 :: Only",
    ],
    python_requires=">=3.5",
)
