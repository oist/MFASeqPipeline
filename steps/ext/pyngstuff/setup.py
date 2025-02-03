import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyngstuff",
    version="0.0.1",
    author="Florian G. Pflug",
    author_email="fgp@phlo.org",
    description="various algorithms and helpers for processing NGS data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/fgp/pyngstuff",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5',
    install_requires=[
        "sortedcontainers",
        "recordclass",
        "wrapt",
        "regex",
        "pysam",
        "progress"
    ]
)
