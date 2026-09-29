# Background materials and research positioning

The local `papers/` directory contains downloaded PDFs; `sources.json` records
URLs, status, file sizes and SHA-256 hashes. The PDFs are not redistributed in
the public repository. Run `python scripts/gather_background.py` to obtain them.
The metadata repository is downloaded separately and excluded from this repository.

| Work | Why it matters | Difference from the present experiment |
|---|---|---|
| [cryoDRGN, Nature Methods 2021](https://doi.org/10.1038/s41592-020-01049-4) | Latent-conditioned Fourier/Hartley coordinate decoder for heterogeneity | Our main experiment fits one consensus volume and uses explicit Gaussian spectral kernels |
| [cryoDRGN2, ICCV 2021](https://openaccess.thecvf.com/content/ICCV2021/html/Zhong_CryoDRGN2_Ab_Initio_Neural_Reconstruction_of_3D_Protein_Structures_From_ICCV_2021_paper.html) | Joint structure/pose inference through hierarchical search | Supplied poses are used here; results must not be labeled ab initio |
| [E2GMM, 2021](https://doi.org/10.1038/s41592-021-01220-5) | Gaussian-mixture cryo-EM heterogeneity already predates recent graphics splatting | It models real-space density mixtures, not nonzero-frequency complex Gaussian kernels |
| [CryoSPIRE](https://arxiv.org/abs/2506.09063) | Part-based hierarchical Gaussian heterogeneity model | Part discovery and heterogeneity are outside our measured experiment |
| [CryoGS/CryoSplat](https://arxiv.org/abs/2508.04929) | Real-space anisotropic splatting adapted to cryo-EM; current arXiv PDF is the ICLR 2026 CryoSplat version | We use an analytic Fourier-plane restriction rather than real-space rasterization |
| [GEM](https://arxiv.org/abs/2509.25075) | Efficient real-space Gaussian reconstruction | Gaussian reconstruction itself is not novel; Fourier transforms are not inherently lossy |
| [CryoBench](https://papers.nips.cc/paper_files/paper/2024/file/a2ef5ba272df8f168dc38037cc946be0-Paper-Datasets_and_Benchmarks_Track.pdf) | Heterogeneous reconstruction benchmarks | Our three accessions are real experimental particle datasets, not the simulated CryoBench suite |
| [Scheres and Chen, 2012](https://doi.org/10.1038/nmeth.2115) | Independent half-set refinement and overfitting prevention | Fixed consensus poses make our half-map FSC conditional rather than fully gold standard |
| [Chen et al., 2013](https://doi.org/10.1016/j.ultramic.2013.06.004) | Noise-substitution tests for inflated FSC | Unmasked FSC and randomized-noise controls are used as transparent checks |

## Central mathematical distinction

A displaced **spatial** Gaussian has a zero-centered Fourier Gaussian envelope
multiplied by a phase ramp. A displaced **spectral** Gaussian has a modulated
Gaussian inverse transform. The proposed conjugate pair is real in spatial
coordinates but can have negative lobes; it is not a positive density mixture.
Restricting a spectral Gaussian to a central plane requires precision restriction
and a plane-distance attenuation factor. Projecting its covariance as if it were
a real-space line integral gives the wrong operation.

The measured stationary lattice case is Gaussian radial-basis regression.
Adaptive centers, anisotropic covariance, tile culling, ab initio pose search,
and heterogeneity cannot be claimed as validated unless separately executed.
No exhaustive novelty claim is made.
