from setuptools import setup, find_packages

with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name             = "KOMETA",
    version          = "1.0.0",
    author           = "Dr. Mohsen Doust Mohammadi",
    author_email     = "",
    description      = "KMC Organic Molecular Exciton Transport Analysis",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url              = "https://github.com/megamsndm/KOMETA",
    packages         = find_packages(),
    python_requires  = ">=3.10",
    install_requires = [
        "numpy>=1.20",
        "scipy>=1.6",
        "matplotlib>=3.3",
        "PyQt5>=5.12",
        "h5py>=3.0",
        "scikit-learn>=0.24",
    ],
    entry_points = {
        "console_scripts": [
            "kometa=kometa_gui:main",
        ],
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Intended Audience :: Science/Research",
    ],
    keywords = [
        "KMC", "kinetic Monte Carlo", "organic semiconductors",
        "charge transport", "triplet diffusion", "TTA upconversion",
        "OLED", "exciton", "Marcus rate", "organic electronics",
        "PEGASUS", "Miller-Abrahams", "Forster", "Holstein",
    ],
)
