import numpy as np
from fourier_splats.uq_pose_cone_probe import InterpolatedCellFourier
from fourier_splats.uq_continuous import cell_forward


def exact(k,coefficients,box):
    values=cell_forward(k[None],np.ones((1,len(k))),coefficients,box,1.).reshape(1,-1)[0]
    return values[:len(k)]+1j*values[len(k):]


def test_padded_fft_cell_guide_is_exact_on_its_frequency_nodes():
    rng=np.random.default_rng(119);box=8;coefficients=rng.normal(size=box**3)
    k=rng.integers(-8,9,size=(20,3))/4
    guide=InterpolatedCellFourier(coefficients,box,32)
    np.testing.assert_allclose(guide.values(k),exact(k,coefficients,box),rtol=1e-9,atol=1e-10)


def test_interpolation_refinement_and_explicit_noncertification():
    rng=np.random.default_rng(121);box=8;coefficients=rng.normal(size=box**3);k=rng.uniform(-2,2,(20,3))
    truth=exact(k,coefficients,box)
    coarse=InterpolatedCellFourier(coefficients,box,32).values(k)
    fine=InterpolatedCellFourier(coefficients,box,64).values(k)
    assert np.linalg.norm(fine-truth)<np.linalg.norm(coarse-truth)/5
    assert np.linalg.norm(fine-truth)/np.linalg.norm(truth)<1e-4
