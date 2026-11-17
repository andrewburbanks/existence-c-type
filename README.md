# existence-c-type

Supporting code for the paper *Existence of the C-type Renormalisation Two-Cycle*.

## Abstract

> We prove the existence of the C-type renormalisation two-cycle, helping to establish the universality of the C-type route to chaos in families of non-invertible maps of the plane.
>
> Families of two-dimensional non-invertible maps with at least two parameters and critical points arising from a fold singularity (one of the fundamental types, alongside the cusp and the projection of the Whitney umbrella), exhibit a distinct type of critical scaling: the C-type. An accumulation of parameter values leads to an infinite collection of coexisting attracting cycles of periods $4^n$ or $2\cdot 4^n$ for which, asymptotically, quadrupling the period is accompanied by parameter-space and state-space scaling governed by particular universal constants. Kuznetsov et al. provided an explanation for this phenomenon in terms of a stationary orbit of period two of the Renormalisation Group (RG) transformation for period-doubling.
>
> We prove the existence of the renormalisation two-cycle in a Banach space of analytic maps and gain rigorous bounds on the corresponding universal state space scaling constants. This result provides a further step in proving a series of outstanding conjectures concerning distinct universality classes for period-doubling. It extends recent results for unidirectionally-coupled maps (the FS-type) to bidirectionally-coupled maps, and generalises the framework from fixed points to periodic orbits of the corresponding renormalisation operators. It provides a further step in establishing the conjectured picture that the C-type universality class is born from the FS-type class via a period-doubling bifurcation in the dynamics of the RG transformation itself. The proof relies on rigorous (validated) computations to establish that a variant of Newton's method for the two-cycle is a contraction map.
>
> The C-type scaling regularity occurs in a number of dynamical systems of interest, perhaps most notably in biologically-plausible models of nephron blood pressure autoregulation.

## Getting Started

The included `Jupyter` notebook `ctype_proof_n24.ipynb` provides the code and results from running the proof with truncation (homogeneous) degree $N=24$ using crude (fast) outward rounding of basic operations to next and previous `Float64` floats.

### Distributed

The notebook uses the `Distributed` module for multiprocessing. The number of processors to use should be set near the top of the notebook by editing the line that reads, e.g.,

```
addprocs(12)
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

### IEEE rounding modes via LLVM

A more sophisticated (tighter bounds, slower) version also uses `Float64` but makes `LLVM` calls for the relevant `IEEE` directed rounding modes. Note that this will only work with particular versions of `Julia`, typically those before `1.10`, because access to rounding modes has been deprecated in some later versions. This mode is selected via

```
const RR = Float64
const use_crude_rounding = false
```

### IEEE rounding modes and multiprecision via `BigFloat`

Finally, for multiprecision support, it is possible to use `BigFloat`. The relevant mode is selected via

```
const RR = BigFloat; setprecision(BigFloat, 256) #example
```

When using datatypes other than `Float64`, the `use_crude_rounding` flag is ignored.

### Requirements

We recommend at least `Julia 1.10.1` if using crude rounding.

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

## Authors

* Zainab Rahman
* Maria Pickett
* Andrew Burbanks - [andrew.burbanks@port.ac.uk](mailto:andrew.burbanks@port.ac.uk)

## License

This code is licensed under the GNU General Public License version 3 (GPLv3).

## Acknowledgements

Work related to this code was supported by the following mechanisms:
* Andrew Burbanks would like to thank the Isaac Newton Institute for Mathematical Sciences (INI), Cambridge, for support and hospitality during an INI Retreat. This work was supported by EPSRC grant EP/Z000580/1.
* The authors would also like to thank the International Centre for Mathematical Sciences (ICMS), Edinburgh, for support via Workshop Grant 2416-BUR, the London Mathematical Society (LMS) for support via the Scheme~1 Conference Grant~12502, and the Glasgow Mathematical Journal Trust (GMJT) for support via the GMJT ICMS Block Grant 2026/2027.
* Zainab Rahman was supported by a PhD bursary from the School of Mathematics and Physics at the University of Portsmouth, UK, through the REF2021 Quality-related Research (QR) funding allocation.
