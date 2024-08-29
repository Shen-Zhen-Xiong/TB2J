import os


def write_umpdss(cls, path="TB2J_results/UMPDSS/input"):
    """
    Write UMPDSS all input
    include INPUT and HAMILTONIAN_PARAMETERS
    """
    if not os.path.exists(path):
        os.makedirs(path)
    write_umpdss_all(cls, path)


def write_umpdss_input(fname, nspins, traverse):
    """
    Write umpdss INPUT file
    which need fname, nspins and traverse as input
    """
    traverse_lattice = ''
    for i in range(3):
        for j in range(2):
            traverse_lattice += f'{traverse[i][j]}  '
    traverse_lattice = traverse_lattice[:-2]
    with open(fname, "w") as myfile:
        myfile.write(f' seed_zero                    -2340\n')
        myfile.write(f' seed_cauchy                  -4321\n')
        myfile.write(f' seed_weight                  -5000\n')
        myfile.write(f' recursion_weight_scale       0.7d0\n')
        myfile.write(f' if_output_debug              .false.\n')
        myfile.write(f' if_perform_exchange          .true.\n')
        myfile.write(f' energy_log_level             2\n')
        myfile.write(f' para_unit_level              2\n')
        myfile.write(f' adjust_spin_scale_freq       5000\n')
        myfile.write(f' max_probability              0.5d0\n')
        myfile.write(f' min_probability              0.4d0\n')
        myfile.write(f' adjust_beta0                 2.0d0\n')
        myfile.write(f' adjust_beta1                 2.0d0\n')

        myfile.write(f' n_scheme                     2\n')
        myfile.write(f' n_recursion                  10\n')
        myfile.write(f' stop_recursion               9\n')
        myfile.write(f' n_sweep_recursion            100\n')
        myfile.write(f' n_exchange_recursion         100\n')
        myfile.write(f' n_sweep_sample               500\n')
        myfile.write(f' n_exchange_sample            200\n')
        myfile.write(f' n_equil_init                 50000\n')
        myfile.write(f' n_equil_recursion            10000\n')
        myfile.write(f' beta_min                     1.159d-2\n')
        myfile.write(f' beta_max                     1.159d1\n')
        myfile.write(f' if_read_beta                 .false.\n')
        myfile.write(f' if_read_betadis              .false.\n')
        myfile.write(f' if_restart                   .false.\n')
        myfile.write(f' if_read_init_configs         .false.\n')

        myfile.write(f' if_normalize_spin            .true.\n')
        myfile.write(f' if_update_spin               .true.\n')
        myfile.write(f' if_fixed_mag_moment          .true.\n')
        myfile.write(f' spin_dimension               3\n')
        myfile.write(f' spin_restrain                1.00\n')
        myfile.write(f' spin_center                  0.0 0.0 0.00\n')
        myfile.write(f' spin_scale                   1.0d0\n')
        myfile.write(f' para_lattice_number          1\n')
        myfile.write(f' spins_number                 {nspins}\n')

        myfile.write(f' spin_lattice                 10       10        10\n')
        myfile.write(f' traverse_lattice             {traverse_lattice}\n')
        myfile.write(f' if_nospin                    .false.\n')
        myfile.write(f' if_exchange                  .true.\n')
        myfile.write(f' if_anisotropy                .false.\n')
        myfile.write(f' if_landau                    .false.\n')
        myfile.write(f' if_biquadratic               .false.\n')
        myfile.write(f' landau_order                 4\n')


def get_umpdss_parameters(cls):
    """
    Get umpdss exchange interactions
    where only colinear Js are including in umpdss_exchange_dict.
    """
    nspins = sum([1 if i > -1 else 0 for i in cls.index_spin])
    umpdss_exchange_dict = {}
    default_data = 0.000
    for si in range(nspins):
        umpdss_exchange_dict[si] = {}
        for sj in range(nspins):
            umpdss_exchange_dict[si][sj] = {}
            for ja in range(-3, 4, 1):
                umpdss_exchange_dict[si][sj][ja] = {}
                for jb in range(-3, 4, 1):
                    umpdss_exchange_dict[si][sj][ja][jb] = {}
                    for jc in range(-3, 4, 1):
                        umpdss_exchange_dict[si][sj][ja][jb][jc] = default_data

    for key in cls.exchange_Jdict:
        R, ispin, jspin = key
        Jiso = -2 * cls.exchange_Jdict[key] * 1e3
        umpdss_exchange_dict[ispin][jspin][R[0]][R[1]][R[2]] = Jiso
    return umpdss_exchange_dict


def write_umpdss_all(cls, fname_head):
    """
    Write umpdss INPUT and HAMILTONIAN_PARAMETERS for five cases
    for 3D bulk, lattice 3(2,1) are avaiable,
    for 2D materials, lattice 3(2,1) are avaiable.
    """

    umpdss_exchange_dict = get_umpdss_parameters(cls)
    fname_tails = ['_3D111', '_3D222', '_3D333', '_2D220', '_2D330']

    for name_tail in fname_tails:
        fname = os.path.join(fname_head, 'HAMILTONIAN_PARAMETERS' + name_tail)
        with open(fname, "w") as myfile:
            myfile.write(" ****START****\n")
            myfile.write(" #energy_nospin(meV)\n")
            myfile.write("          0.000\n")
            myfile.write(" #Exchange(jc,jb,ja,sj,si)(meV)\n")
            a_range = [-1, 1]
            b_range = [-1, 1]
            c_range = [-1, 1]
            if name_tail == '_3D222':
                a_range = [-2, 2]
                b_range = [-2, 2]
                c_range = [-2, 2]
            elif name_tail == '_3D333':
                a_range = [-3, 3]
                b_range = [-3, 3]
                c_range = [-3, 3]
            elif name_tail == '_2D220':
                a_range = [-2, 2]
                b_range = [-2, 2]
                c_range = [0, 0]
            elif name_tail == '_2D330':
                a_range = [-3, 3]
                b_range = [-3, 3]
                c_range = [0, 0]
            nspins = sum([1 if i > -1 else 0 for i in cls.index_spin])
            traverse = [a_range, b_range, c_range]
            input_name = os.path.join(fname_head, 'INPUT' + name_tail)
            write_umpdss_input(input_name, nspins, traverse)
            for si in range(nspins):
                for sj in range(nspins):
                    for ja in range(a_range[0], a_range[1]+1, 1):
                        for jb in range(b_range[0], b_range[1]+1, 1):
                            for jc in range(c_range[0], c_range[1]+1, 1):
                                flag = (ja == 0 and 
                                        jb == 0 and 
                                        jc == 0 and
                                        si == sj)
                                if (not flag):
                                    Jiso = umpdss_exchange_dict[si][sj][ja][jb][jc]
                                    data = " ".join("{:>9.3f}".format(0.000) 
                                                    for _ in range(8)) + " {:>9.3f}".format(Jiso)
                                    myfile.write(f'{data}\n')

            myfile.write(" #anisotropy(0:2,0:2)(meV)\n")
            output_anisotropy = '    0.000' * 9
            myfile.write(f"{output_anisotropy}\n")
            myfile.write(" #landau_coef(0:3)(meV)\n")
            output_landau = '    0.000' * 4
            myfile.write(f"{output_landau}\n")
            myfile.write(" #biquadratic(jc,jb,ja,sj,si)(meV)\n")
            for si in range(nspins):
                for sj in range(nspins):
                    for ja in range(a_range[0], a_range[1]+1, 1):
                        for jb in range(b_range[0], b_range[1]+1, 1):
                            for jc in range(c_range[0], c_range[1]+1, 1):
                                flag = (ja == 0 and 
                                        jb == 0 and 
                                        jc == 0 and
                                        si == sj)
                                if (not flag):
                                    myfile.write('    0.000\n')
            myfile.write(f'  ****END****\n')