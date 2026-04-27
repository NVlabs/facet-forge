# FacetForge

FacetForge is a reference C++ codebase that provides benchmark implementations of microfacet BSDFs with unbiased multiple scattering.  Unique capabilities:
- Quickly add new NDFs (including full-sphere NDFs) by implementing one simple function (the NDF itself)
- Visible-distribution-of-normals sampling for Student-T NDFs
- Support for biscale roughness: put any BSDF onto the microsurface, including other microsurface BSDFs

The initial release of this codebase was part of:
> __Student-T and Beyond: Practical Tools for Multiple-Scattering BSDFs with General NDFs__  
> [Eugene d'Eon](http://eugenedeon.com)  
> _ACM (__SIGGRAPH__) Talks, August 2023_  
> __[Paper](assets/deon2023student.pdf)&nbsp;/ [Supplemental](assets/deon2023studentsupp.pdf)__

We will update the codebase over time.  Pull requests are encouraged.  For guidelines, see [CONTRIBUTING](CONTRIBUTING.md).

## Usage

### Microsurface

Different from other approaches, in FacetForge a Microsurface is an operator that takes as input a BSDF and an NDF and outputs a new BSDF.  Therefore, there are no RoughDielectric or RoughConductor BSDFs.  To create a rough conductor, for example, you first create a smooth ConductorBRDF, plug it into your chosen NDF, and plug that into the Microsurface class to create a new BRDF:
```
ConductorBRDF micro_brdf(eta, k); // take a smooth conductor BRDF
GGXNDF ndf(&micro_brdf, rough_x, rough_y); // and assign it to microfacets with a GGX distribution
Microsurface macro_brdf(&ndf); // and feed this to a Microsurface operator to create a rough conductor BRDF
```
More examples of rough BSDFs are included in the `test` folder.

### NDFs

NDFs can be added to FacetForge in two ways:
- heightfield NDFs can derive from the `ShapeInvariantNDF` and must implement the `P22` slope distribution (which defines the NDF), sampling of the visible distribution of slopes when both roughnesses are equal to unity, and the cross section as a function of direction over the full sphere
- general full-sphere NDFs can derive from the `NullNDF` and implement the NDF `D` together with a majorant

## Limitations

The primary purpose of the codebase is to implement flexible microfacet BSDFs with general NDFs.  Achieving this goal comes with some limitations (some of which are straightforward to remove, some not), including:
- `pdf()` is not implemented
- NullNDFs will suffer crippling inefficiency for very low roughness (analogous to null scattering through a mostly empty inhomogeneous medium with a very large majorant)
- Analytic `eval` for single-scattering and specular facets is not currently implemented
- Polarization is not currently supported
- There is no notion of spectrum or color - radiance is monochromatic `double`

## Assumptions

The code has three parts:
- C++ implementation (generally portable) in the `include` folder
  - assumes `drand48()`
  - tested on Mac OS Arm M1 with `clang++`
  - assumes `std::mt19937` for gamma random variates (for Student-T NDF sampling)
- Mathematica tests (described below) in the `test` folder
- Python test and plotting front ends next to the C++ tests in the `test` folder

## Running the Tests

The `test` folder contains C++ test drivers, Mathematica notebooks, and Python plotting front ends that compare `eval` and `sample` for various BSDFs and NDFs.

The Python tests are intended to be the simplest reproducible workflow.  Each Python file sits next to the matching `.cpp` and `.nb` file, compiles the C++ driver, runs it with defaults matching the notebook parameters, and writes a `.png` comparison plot.  The 2D comparison plots include gamma-corrected heatmaps plus theta and phi marginals.  The 1D sigma tests compare discrete/Monte Carlo outputs to deterministic or analytic reference curves where available.

The Python plotting scripts use common scientific Python packages:
```
python3 -m pip install numpy matplotlib scipy
```

For example:
```
python3 test/lambert/test_lambert_eval_sample.py
python3 test/NDFs/test_ST_sigma.py
```

To run the Python test suite and build a scrollable HTML image report:
```
python3 test/run_all_python_tests.py --clean
```
The report is written to `test/_python_report/latest/index.html`.  The top of the page links to every test with blue pass chips and red failure chips; each section embeds the generated `.png` and the numeric summary printed by the wrapper.  Use `--only '*dielectric*'` to run a subset or `--jobs 4` to run independent tests concurrently.

The notebook-matched default sample counts are intentionally large.  For quick checks, override them on the command line:
```
python3 test/lambert/test_lambert_eval_sample.py --samples 100000 --eval-samples 4
python3 test/NDFs/test_ST_sigma.py --du 0.05
```

Use `--compiler clang++` or `--compiler g++` to select the C++ compiler, and `--output path/to/plot.png` to choose the output image path.

The original Mathematica notebooks remain useful for interactive exploration.  To use the notebooks:
- In Mathematica `SetDirectory[]` to the `facet-forge` root dir where you cloned to
- (optional): edit the arguments of the test to vary roughness parameters, iors, incidence angle, numbers of samples
- run the full notebook, which will:
  - compile the required C++ test using `clang++`
  - execute the test, dumping the output to a `.txt` file
  - load the output and compare 2D histograms and 1D curves that compare `eval` and `sample`

## Thanks

This codebase is an extension of Eric Heitz's original implementation from the 2016 SIGGRAPH paper that introduced multiple scattering to microfacet theory:
https://eheitzresearch.wordpress.com/240-2/

## License and Citation

```bibtex
@incollection{deon2023,
    author = {Eugene d'Eon},
    title = {Student-T and Beyond: Practical Tools for Multiple-Scattering BSDFs with General NDFs},
    booktitle={ACM SIGGRAPH 2023 Talks},
    pages={1--2},
    year={2023},
    url={https://doi.org/10.1145/3587421.3595417}
}
```

Copyright © 2023, NVIDIA Corporation. All rights reserved.

This code is made available under the Apache-2.0 license.
