"""Differentiable anisotropic paired Fourier kernels (dense reference renderer).

This reference validates the general proposed representation. The main real-data
benchmark uses the sparse fixed-geometry solver in basis.py, not this renderer.
"""

import torch
from torch import nn


class FourierGaussianPairs(nn.Module):
    def __init__(self, centers, coefficients, precision_cholesky):
        super().__init__()
        self.centers = nn.Parameter(torch.as_tensor(centers).clone())
        coefficients = torch.as_tensor(coefficients)
        self.coefficients = nn.Parameter(
            torch.stack([coefficients.real, coefficients.imag], -1).clone()
        )
        L = torch.as_tensor(precision_cholesky).clone()
        diagonal = torch.diagonal(L, dim1=-2, dim2=-1)
        self.log_diagonal = nn.Parameter(torch.log(diagonal))
        self.lower = nn.Parameter(torch.stack([L[:, 1, 0], L[:, 2, 0], L[:, 2, 1]], -1))

    def precision_factor(self):
        L = torch.diag_embed(torch.exp(self.log_diagonal))
        L[:, 1, 0] = self.lower[:, 0]
        L[:, 2, 0] = self.lower[:, 1]
        L[:, 2, 1] = self.lower[:, 2]
        return L

    def forward(self, k):
        """Return real/imag channels; no complex MPS operations are required."""
        L = self.precision_factor()
        dp = k[:, None, :] - self.centers[None, :, :]
        dm = k[:, None, :] + self.centers[None, :, :]
        gp = torch.exp(-0.5 * torch.einsum("nki,kij->nkj", dp, L).square().sum(-1))
        gm = torch.exp(-0.5 * torch.einsum("nki,kij->nkj", dm, L).square().sum(-1))
        real = (gp + gm) @ self.coefficients[:, 0]
        imag = (gp - gm) @ self.coefficients[:, 1]
        return torch.stack([real, imag], -1)
