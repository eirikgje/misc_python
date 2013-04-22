from __future__ import division
import numpy as np
from scipy.special import gamma

def normalize_1d(dist, dx):
    return dist / (np.sum(dist) * dx)

def inverse_gamma(x, alpha, beta):
    dx = x[1] - x[0]
    dist = -(alpha + 1) * np.log(x) - beta / x
    dist = np.exp(dist - np.max(dist))
    return dist / np.sum(dist) / dx

def gaussian(x, mean, std):
    dx = x[1] - x[0]
    dist = -np.log(std * np.sqrt(2 * np.pi)) - 0.5 * ((x - mean) / std) ** 2
    dist = np.exp(dist - np.max(dist))
    return dist / np.sum(dist) / dx
