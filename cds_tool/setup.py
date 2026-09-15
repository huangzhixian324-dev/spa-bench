from setuptools import setup

setup(
    name="circularity-detection-score",
    version="1.2.1",
    description="Detect circular dependencies between benchmark endpoints and prediction features in transcriptomic studies",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="SPATBench Team",
    url="https://github.com/spatbench/cds",
    py_modules=["cds"],
    install_requires=["numpy>=1.21", "scipy>=1.7", "scikit-learn>=1.0", "pandas>=1.3"],
    entry_points={
        "console_scripts": [
            "cds=cds:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
)
