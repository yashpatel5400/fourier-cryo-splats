# Alignment-moment follow-up after the review-3 snapshot

These targeted readings form a follow-up queue; they are not silently added to the 131-candidate reviewed manuscript. The running review and its evidence remain unchanged.

[Baldwin and Penczek (2005)](https://pubmed.ncbi.nlm.nih.gov/15866744/) already estimate translational, rotational and Gaussian-noise variances through Fourier-harmonic image averages and intensities with an unknown pattern. Only the author abstract was checked; it does not establish a usable three-dimensional experimental pose-confidence ball for this project.

[Park et al. (2011)](https://rpk.lcsr.jhu.edu/wp-content/uploads/2014/08/Classaverage_ijrr_final.pdf) model alignment blur through convolution on SE(2). Their chosen alignment matches image mass centers and principal axes. The rotational derivation uses a small-angle approximation; treating its distribution as Gaussian additionally neglects a product term judged small in their example. Their experimental data are negative-stain receptor images with corrected CTFs, and their noise model includes spatial correlations. Validation compares simulated and experimental class averages with convolution predictions. It is not a calibrated confidence procedure for a three-dimensional density feature. Our inference from these passages is that averaging the nonlinear pose operator or modeling correlated background noise is established prior work; neither alone is a new contribution. Transfer to estimated three-dimensional slice geometry and post-alignment inference would require a separate argument.

The source record states exactly which passages were read and which downloads failed. No new experiment or claimed uncertainty improvement follows from this note.
