import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="asali",
    version="1.0.0",
    author="StefanoRebughini",
    author_email="ste.rebu@gmail.com",
    description="ASALI: Modeling and beyond",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="www.catalyticfoam.polimi.it/asali",
    include_package_data=True,
    packages=["asali",
              "asali.gasproperties",
              "asali.reactors",
              "asali.unittest",
              "asali.unittest.BatchReactor",
              "asali.unittest.CstrReactor",
              "asali.unittest.Heterogeneous1DReactor",
              "asali.unittest.PseudoHomogeneous1DReactor",
              "asali.unittest.utils",
              "asali.database"],
    classifiers=["Programming Language :: Python :: 3", "Operating System :: OS Independent"],
    install_requires=[],
    python_requires='>=3.6'
)
