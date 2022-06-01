#!/usr/bin/env python

'''
Hartree-Fock/DFT with k-points sampling for all-electron calculations

GDF (Gaussian density fitting), MDF (mixed density fitting), RSGDF
(range-separated Gaussian density fitting), or RS-JK builder
can be used in all electron calculations. They are more efficient than the
default SCF JK builder.
'''

import numpy
from pyscf.pbc import gto, scf, dft

cell = gto.M(
    a = '''
        2.8664000034         0.0000000000         0.0000000000
        0.0000000000         2.8664000034         0.0000000000
        0.0000000000         0.0000000000         2.8664000034
''',
    atom='''
Fe 0.000000000000   0.000000000000   0.000000000000
Fe 1.433200002      1.433200002  1.433200002''',
    basis = 'gth-szv-molopt-sr',
    pseudo= 'gth-pbe',
    spin = 8,
    verbose = 4,
)


#
# Mixed density fitting
#
#kmf = scf.KRHF(cell, kpts).mix_density_fit()
# In the MDF scheme, modifying the default mesh for PWs to reduce the cost
# The default mesh for PWs is a very dense-grid scheme which is automatically
# generated based on the AO basis. It is often not necessary to use dense grid
# for MDF method.
#kmf.with_df.mesh = [10,10,10]
#kmf.kernel()
#
for nk in range(2,6):
    nk = [nk,nk,nk]  # 4 k-poins for each axis, 4^3=64 kpts in total
    kpts = cell.make_kpts(nk)
    print(kpts)
    kmf = dft.KUKS(cell, kpts).rs_density_fit()
    kmf = scf.addons.smearing_(kmf, sigma=0.1, method='gauss')
    kmf.xc = 'pbe'
    kmf.kernel()

exit()
#
# RS-JK builder is efficient for large number of k-points
#
kmf = scf.KRHF(cell, kpts).jk_method('RS')
kmf.kernel()

#
# Second order SCF solver can be used in the PBC SCF code the same way in the
# molecular calculation.  Note second order SCF algorithm does not support
# smearing method.
#
mf = scf.KRHF(cell, kpts).density_fit()
mf = mf.newton()
mf.kernel()
