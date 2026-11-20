# existence-c-type

Supporting code for the paper *Existence of the C-type renormalisation two-cycle*. Provides rigorous enclosures for operations on balls of functions analytic on bidiscs. Applies this to prove the existence of the c-type two-cycle of the renormalisation operator for period-doubling in pairs of maps of the plane.

## Authors

* Zainab Rahman [ORCiD:0009-0000-6028-8358](https://orcid.org/0009-0000-6028-8358)
* Maria Pickett [ORCiD:0000-0002-3601-1989](https://orcid.org/0000-0002-3601-1989)
* Andrew Burbanks [ORCiD:0000-0003-0685-6670](https://orcid.org/0000-0003-0685-6670) (*)

(*) Corresponding Author - [andrew.burbanks@port.ac.uk](mailto:andrew.burbanks@port.ac.uk)

## Abstract

We prove the existence of the C-type renormalisation two-cycle, helping to establish the universality of the C-type route to chaos in families of non-invertible maps of the plane.

Families of two-dimensional non-invertible maps with at least two parameters and critical points arising from a fold singularity (one of the fundamental types, alongside the cusp and the projection of the Whitney umbrella), exhibit a distinct type of critical scaling: the C-type. An accumulation of parameter values leads to an infinite collection of coexisting attracting cycles of periods $4^n$ or $2\cdot 4^n$ for which, asymptotically, quadrupling the period is accompanied by parameter-space and state-space scaling governed by particular universal constants. Kuznetsov et al. provided an explanation for this phenomenon in terms of a stationary orbit of period two of the Renormalisation Group (RG) transformation for period-doubling.

We prove the existence of the renormalisation two-cycle in a Banach space of analytic maps and gain rigorous bounds on the corresponding universal state space scaling constants. This result provides a further step in proving a series of outstanding conjectures concerning distinct universality classes for period-doubling. It extends recent results for unidirectionally-coupled maps (the FS-type) to bidirectionally-coupled maps, and generalises the framework from fixed points to periodic orbits of the corresponding renormalisation operators. It provides a further step in establishing the conjectured picture that the C-type universality class is born from the FS-type class via a period-doubling bifurcation in the dynamics of the RG transformation itself. The proof relies on rigorous (validated) computations to establish that a variant of Newton's method for the two-cycle is a contraction map.

The C-type scaling regularity occurs in a number of dynamical systems of interest, perhaps most notably in biologically-plausible models of nephron blood pressure autoregulation.

## Manifest

* `dat/` - data files produced by successful run.
* `fig/` - figures produced for publication.
* `rigorous/` - Julia code adapting the framework of Eckmann, Koch, and Wittwer. This provides rigorous bounds for operations on balls of functions analytic on bidiscs.
* `ctype_proof_n24.ipynb` - Jupyter notebook that reproduces the proof.
* `ctype_proof_n24.html` - static web page version of the notebook; download and open in a browser.

## Requirements

We have used `Julia 1.10.1` for the results included here. We have also tested on the latest stable version `Julia 1.13.4`. The following `Julia` packages are used:
* `Distributed`, `Serialization`
* `Plots`, `LaTeXStrings`
* `GenericSchur`
* `IJulia` (to ensure the relevant kernel is available via Jupyter)

## Getting Started

The included `Jupyter` notebook `ctype_proof_n24.ipynb` provides the code and results from running the proof with (homogeneous) truncation degree `N=24` using crude (fast) outward rounding of upper and lower interval bounds on basic operations by 1ULP to the next and previous `Float64` floats.

Note: the file `ctype_proof_n24.html` is a static web page version of the notebook (and its output) that can be downloaded and opened with a web browser, so that Jupyter (and Julia) do not need to be installed in order to see the results.

### Distributed

The notebook uses the `Distributed` module for multiprocessing. The number of processors to use should be set near the top of the notebook by editing the line that reads, e.g.,

```
addprocs(12)
```

### Truncation degree

The homogeneous truncation degree for the polynomial part of the space is set using, for example

```
const N::Int64 = 24
```

### Crude rounding

The rounding model enabled in the demonstration notebook works with a wide variety of `Julia` versions. This mode is selected using the following flags near the top of the notebook (here, `RR` denotes the datatype to use for Representable Reals):

```
const RR = Float64
const use_crude_rounding = true
```

### Using previously-computed results

Data files stored during this run are provided in the `dat/` folder. Loading these files, rather than repeating computations, is achieved using

```
const use_precomputed = true
```

Note: precomputed results are provided for `N=24` using `RR = Float64` and `use_crude_rounding = true`.

### Reproducing figures

Some of the figures in the paper are produced directly using the notebook. These are written to the `fig/` folder. The following flag is used to enable reproduction of the figures

```
const draw_figures = true
```

### IEEE rounding modes and multiprecision via `BigFloat`

For multiprecision support, and/or to use directed rounding modes, it is possible to use `BigFloat`. The relevant mode is selected via

```
const RR = BigFloat; setprecision(BigFloat, 64) #example
```

Note that when using datatypes other than `Float64`, the `use_crude_rounding` flag is ignored.

### IEEE rounding modes via LLVM

Code is also provided that uses `Float64` but makes `LLVM` calls for the relevant `IEEE` directed rounding modes. Note that this will only work with particular versions of `Julia`, typically those before `1.10`, because access to rounding modes has been deprecated in some later versions. This mode is selected via

```
const RR = Float64
const use_crude_rounding = false
```

## Publications

The following preprint made use of this code.

__Existence of the C-type renormalisation two-cycle__<br />
Zainab Rahman, Maria Pickett, Andrew Burbanks<br />
2026. arXiv 2606.06274, math.DS<br />
DOI: [https://doi.org/10.48550/arXiv.2606.06274](https://doi.org/10.48550/arXiv.2606.06274)<br />
URL: [https://arxiv.org/abs/2606.06274](https://arxiv.org/abs/2606.06274) 

## Contributing

Ideas and comments are welcome. Please email the corresponding author.

## Notes

As indicated above, when using `Float64` as the type for Representable Reals, LLVM access to rounding modes will only work in versions of Julia that provide access to it; Set `use_crude_rounding = true` in order to avoid LLVM use if needed. This does not affect `BigFloat` which provides its own directed rounding modes.

## Change log

* 2026-11-18 Fresh repo with readme.md files, data, figures, .html version.
* 2026-06-04 Initial preprint on arXiv.
* 2025-10-01 Rigorous framework and proof code initial version.

## License

This code is licensed under the GNU General Public License version 3 (GPLv3).

## Acknowledgements

Work related to this code was supported by the following mechanisms:
* Andrew Burbanks would like to thank the Isaac Newton Institute for Mathematical Sciences (INI), Cambridge, for support and hospitality during an INI Retreat. This work was supported by EPSRC grant EP/Z000580/1.
* The authors would also like to thank the International Centre for Mathematical Sciences (ICMS), Edinburgh, for support via Workshop Grant 2416-BUR, the London Mathematical Society (LMS) for support via the Scheme 1 Conference Grant 12502, and the Glasgow Mathematical Journal Trust (GMJT) for support via the GMJT ICMS Block Grant 2026/2027.
* Zainab Rahman was supported by a PhD bursary from the School of Mathematics and Physics at the University of Portsmouth, UK, through the REF2021 Quality-related Research (QR) funding allocation.

## References

The following references were useful in creating the code:
* _A Computer-Assisted Proof of Universality for Area-Preserving Maps._ Jean Pierre Eckmann, Hans Koch, and Peter Wittwer. Memoirs of the American Mathematical Society, volume 47, number 289, ISBN 1-4704-0699-3; Providence, R.I., USA : American Mathematical Society 1984.
* _Existence of the FS-type renormalisation ﬁxed point for unidirectionally-coupled pairs of maps._ Burbanks, A. D. & Osbaldestin, A. H., 20 Apr 2023, In: Journal of Physics A: Mathematical and Theoretical. 56, 17 p., 195202. DOI: [10.1088/1751-8121/acc7dd](https://doi.org/10.1088/1751-8121/acc7dd)
* _andrewburbanks/existence-fs: Existence of the FS-type renormalisation fixed point [Software]_ Burbanks, A. (Creator), Zenodo, 3 Oct 2022. DOI: [10.5281/zenodo.7139006](https://doi.org/10.5281/zenodo.7139006)
* _Rigorous computer-assisted bounds on the period doubling renormalisation fixed point and eigenfunctions in maps with critical point of degree 4._ Burbanks, A., Osbaldestin, A. & Thurlby, J., 30 Nov 2021, In: Journal of Mathematical Physics. 62, 11, 18 p., 112701. DOI: [10.1063/5.0054823](https://doi.org/10.1063/5.0054823)
* _Dataset for "Rigorous computer-assisted bounds on the period doubling renormalisation fixed point and eigenfunctions in maps with critical point of degree 4"._ Burbanks, A., Zenodo, 28 Oct 2021, DOI: [10.5281/zenodo.5608449](https://doi.org/10.5281/zenodo.5608449)

