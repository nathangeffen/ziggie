import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="epi",
    version="0.0.1",
    author="Nathan Geffen",
    author_email="nathangeffen@gmail.com",
    description="Infectious disease modelling package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nathangeffen/epi",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
    ],
    python_requires='>=3.5',
)
