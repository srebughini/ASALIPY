import setuptools

setuptools.setup(
    name="asali",
    version="0.1",
    author="ASALI",
    author_email="ste.rebu@outlook.it",
    description="Library to model chemical reactors",
    url="https://github.com/srebughini/ASALIPY",
    classifiers=["Programming Language :: Python :: 3", "Operating System :: OS Independent"],
    package_dir={"asali": "src"},
    packages=["asali.utils", "asali.plotters", "asali.reactors"],
    include_package_data=True
)