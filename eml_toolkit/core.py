import mpmath as mp

mp.mp.dps = 60

ONE = mp.mpf(1)

def eml(x, y):
    """eml(x, y) = exp(x) - ln(y)"""
    return mp.exp(x) - mp.log(y)

def eml_star(x, y):
    """eml★(x, y) = exp(x) - ln(conj(y))"""
    return mp.exp(x) - mp.log(mp.conj(y))

def fold_to_strip(z):
    """Folding continu Im(z) → [-π, π) sans saut discontinu"""
    re = mp.re(z)
    im = mp.im(z)
    folded_im = mp.fmod(im + mp.pi, mp.mpf(2)*mp.pi) - mp.pi
    return mp.mpc(re, folded_im)

def E_neg(z):
    return eml(eml(ONE, eml(eml(ONE, mp.mpf(0)), ONE)), eml(z, ONE))

def E_sub(a, b):
    return eml(a, eml(b, ONE))

def E_add(a, b):
    return E_sub(a, E_neg(b))

def E_mul(a, b):
    return eml(E_add(eml(a, ONE), eml(b, ONE)), ONE)

def E_conj(z):
    return ONE - eml_star(mp.mpc(0), eml(z, ONE))

def real_part(z):
    return mp.re(z)

def imag_part(z):
    return mp.im(z)

def modulus_squared(z):
    return eml(z, eml_star(z, ONE))

def modulus(z):
    return mp.sqrt(modulus_squared(z))

__all__ = ['eml', 'eml_star', 'fold_to_strip', 'E_add', 'E_mul', 'E_conj',
           'real_part', 'imag_part', 'modulus_squared', 'modulus', 'ONE']
