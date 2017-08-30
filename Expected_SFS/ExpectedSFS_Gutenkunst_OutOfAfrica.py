# November 28, 2016

# Author: Beichman, Annabel
# Contact: annabel.beichman@gmail.com
# Purpose: Python code modified from Gutenkunst et al. (2009), Supplementary Information
# to get expected 3-D joint SFS under the Out of Africa Model
# (this SFS is then marginalized to get the SFSes of each population)

def OutOfAfrica ((nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs, mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs),(n1,n2,n3), pts):
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, TAf, nu= nuAf)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, TB, nu1=nuAf, nu2 =nuB,m12 =mAfB, m21 = mAfB)
    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    nuEu_func = lambda t: nuEu0 *(nuEu / nuEu0)**(t/ TEuAs)
    nuAs_func = lambda t: nuAs0 *(nuAs / nuAs0)**(t/ TEuAs)
    phi = Integration.three_pops(phi, xx, TEuAs, nu1 =nuAf,nu2 = nuEu_func, nu3= nuAs_func,m12 =mAfEu, m13 =mAfAs, m21 =mAfEu,m23 =mEuAs, m31 =mAfAs, m32 = mEuAs)
    fs = Spectrum.from_phi(phi, (n1,n2,n3), (xx,xx,xx))
    return fs

gutenkunst_extrap_function = Numerics.make_extrap_log_func(OutOfAfrica)

# Parameters from Gutenkunst et al. (2009) SI pg 12
# (nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs, mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs)
# Note: dadi model has times as time from present, but should be intervals
gutenkunst_params=(1.68, 0.287, 0.129, 3.74, 0.07, 7.29, 3.65, 0.44, 0.28, 1.4, (0.607-.396), (0.396-.058), 0.058)
