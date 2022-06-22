from __future__ import division, absolute_import, print_function

def calculate():
    import psi4
    import resp

    mol = psi4.geometry("""
-1 1

    C            0.613146075611     1.296817574710    -0.001408652073
    C           -0.054498764942    -0.101727893662    -0.001088564526
    H           -0.135232499598     2.094735817787    -0.000936517573
    H            1.255589748430     1.392531005997    -0.884205440087
    H            1.256637717162     1.392473067087     0.880634402951
    O           -1.324128058114    -0.123476765848     0.001602292852
    O            0.755237071679    -1.080593442882     0.000555231991

    """)

    mol.update_geometry()

    options = {'N_VDW_LAYERS'       : 4,
               'VDW_SCALE_FACTOR'   : 1.4,
               'VDW_INCREMENT'      : 0.2,
               'VDW_POINT_DENSITY'  : 20.0,
               'resp_a'             : 0.0005,
               'RESP_B'             : 0.1,
               'BASIS_ESP':'3-21G',
               'METHOD_ESP':'HF',
               'RADIUS':{'BR':1.97,'I':2.19}
    }

    # Call for first stage fit
    charges1 = resp.resp([mol], [options])

    
calculate()
