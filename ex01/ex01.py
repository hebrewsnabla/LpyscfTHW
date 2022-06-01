#!/usr/bin/env python

'''Fermi-Dirac or Gaussian smearing for PBC SCF calculation'''

import numpy
from pyscf.pbc import gto, scf
from pyscf import lib
lib.num_threads(1)

cell = gto.Cell()
cell.atom = '''
O 0 0 0
'''
cell.basis = 'ccpvdz'
cell.a = numpy.eye(3) * 8
cell.spin=2
cell.verbose = 4
cell.build()

#
# Use scf.addons.smearing_ function to modify the PBC (gamma-point or k-points)
# SCF object
#
#nks = [2,1,1]
mf = scf.KUKS(cell, 
        #cell.make_kpts(nks),
        xc='PBE'
        ).density_fit()
mf = scf.addons.smearing_(mf, 
        sigma=0.01, 
        method='gauss'
        )
mf.kernel()
print('Entropy = %s' % mf.entropy)
print('Energy without entropy = %s  # E' % mf.e_tot)
print('Free energy = %s             # E - sigma*entropy' % mf.e_free)
print('Zero temperature energy = %s # E - 0.5*sigma*entropy = E0' % ((mf.e_tot+mf.e_free)/2))
mo_occ = mf.get_occ(mf.mo_energy, mf.mo_coeff)
print('occ alpha/beta')
print(mo_occ[0])
print(mo_occ[1])

#exit()
#
# The smearing method and parameters can be modified at runtime
#
mf.sigma = .1
mf.method = 'gauss'
mf.max_cycle = 2
mf.kernel()
print('Entropy = %s' % mf.entropy)

mf.sigma = .01
mf.method = 'fermi'
mf.max_cycle = 50
mf.kernel()
print('Entropy = %s' % mf.entropy)

