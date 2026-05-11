from src.yall.state import init_states
from scidatacontainer import Container

from yall import Coupling
from yall.ameli import update, AMELI_PATH

if __name__ == "__main__":
    radial = { "base": 327.39, "H1/2": 68576.05, "H1/4": 49972.76, "H1/6": 32415.29, "H2": 728.18,
        "H3/0": 16.99, "H3/1": -417.98, "H3/2": 1371, "H5fix": 0.19, "H6fix": 1.67 }

    num = 2
    update(num)
    coupling = Coupling.SLJ
    name = "H1/2"

    # with Lanthanide(num, coupling, radial) as ion:
    #     matrix_old = ion.matrix(name, coupling).array
    # matrix_new = read_f_matrix(num, coupling, name)
    # print(matrix_old.shape, matrix_old.dtype)
    # print(matrix_new.shape, matrix_new.dtype)
    # print(f"New matrix {name} is identical: {np.allclose(matrix_old, matrix_new, atol=1e-6)}")
    # np.set_printoptions(formatter=cast(Any, {'float': '{: .3f}'.format}), linewidth=120)
    # #print(matrix_old)
    # #print(matrix_new)

    config_name = f"f{num}"
    dc = Container(file=str(AMELI_PATH / config_name/ "transform.zdc"))

    #t = read_transform(num)
    s = init_states(num)