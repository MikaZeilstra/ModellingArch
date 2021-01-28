import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ModellingArch",
    version="0.0.8",
    author="MikaZ",
    author_email="M.J.Zeilstra@student.tudelft.nl",
    description="A package for modelling the photocycle of arch",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MikaZeilstra/ModellingArch",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    install_requires = [ 
        'matplotlib==3.3.3',
        'sympy==1.7.1',
        'scipy==1.5.4'
    ]
)