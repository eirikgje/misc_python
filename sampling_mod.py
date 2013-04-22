from __future__ import division
import numpy as np

def sample_cls_given_sigma(sigma, l, numsamps):
    sigmat = (2*l + 1) * 2 * np.pi / (l * (l + 1)) * np.array([[sigma[0], sigma[2]], [sigma[2], sigma[1]]])
    inv_sigmat = np.linalg.inv(sigmat)
    inv_sigmat = np.linalg.cholesky(inv_sigmat)
    clmat = np.zeros((2, 2, numsamps))
    for m in range(numsamps):
        for k in range(2 * l - 2):
            eta = np.random.normal(size=2)
            eta = np.dot(inv_sigmat, eta)
            for i in range(2):
                for j in range(2):
                    clmat[i, j, m] += eta[i] * eta[j]
        clmat[:, :, m] = np.linalg.inv(clmat[:, :, m])
    clmat *= l * (l + 1) / (2 * np.pi)
    return clmat
