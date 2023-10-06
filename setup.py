from setuptools import setup

# so mypy can see version
__version__ = "0.0.0"

exec(open("sanclone/version.py").read())

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="sanclone",
    version=__version__,
    description="Molecular cloning agent",
    author="Andrew White",
    author_email="andrew@futurehouse.org",
    url="https://github.com/whitead/sanclone",
    license="MIT",
    packages=["sanclone", "sanclone.tools"],
    install_requires=["langchain", "biopython"],
    test_suite="tests",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)
