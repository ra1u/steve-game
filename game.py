import argparse
from collections import defaultdict

import numpy as np
from ortools.math_opt.io.python import mps_converter
from ortools.math_opt.python import mathopt
from scipy.sparse import dok_array


def gen_cons(n):
    dd = defaultdict(lambda: defaultdict(float))
    for a in range(n):
        for b in range(a, n):
            node_up = (a, b)
            for x in range(a, b + 1):
                node_move = (node_up, x)
                dd[node_up][node_move] += 1.0
                if x > a:
                    # balmer sais - number is less than x
                    node_down = (a, x - 1)
                    dd[node_down][node_move] -= 1.0
                if x < b:
                    # balmer sais - number is more than x
                    node_down = (x + 1, b)
                    dd[node_down][node_move] -= 1.0
    return dd


def node_enums(n):
    ne = [((a, b), x) for a in range(n) for b in range(a, n) for x in range(a, b + 1)]
    nes = set(c_node_symetric(n, nn) for nn in ne)
    nesx = list(nes)
    nesx.sort()
    return nesx


def c_node_symetric(n, node):
    # returns either ((a,b),k) or ((n-b,n-a),n-k)
    ((a, b), k) = node
    asy = n - b - 1
    bsy = n - a - 1
    if a + b < asy + bsy:
        return ((a, b), k)
    else:
        ksy = n - k - 1
        return ((asy, bsy), ksy)


def b_node_symetric(n, node):
    # returns either (node)  or (n - node)
    return min(node, n - node - 1)


def solve(args):
    n = args.count
    c_nodes = node_enums(n)
    c_nodes_i = dict((v, i) for i, v in enumerate(c_nodes))
    cons = gen_cons(n)

    n_bnodes = n // 2 + (n % 2)
    # Build A matrix
    A = dok_array((n_bnodes, len(c_nodes)))
    for a in range(n):
        for b in range(a, n):
            for x in range(a, b + 1):
                node_name = c_node_symetric(n, ((a, b), x)) 
                i = c_nodes_i[node_name]
                for j in range(a, b + 1):
                    js = b_node_symetric(n, j)
                    A[js, i] += 1.0

    # C and c
    c_columns = len(cons)
    C = dok_array((c_columns, len(c_nodes)))
    cc = dok_array((c_columns, 1))
    kk1 = (0, n - 1)
    for i, (k, d) in enumerate(cons.items()):
        if k == kk1:
            cc[i, 0] = 1.0
        for dk, dv in d.items():
            dki = c_node_symetric(n, dk)
            j = c_nodes_i[dki]
            C[i, j] += dv

    # build B
    B = dok_array((1, n_bnodes))
    for j in range(n):
        B[0, b_node_symetric(n, j)] += 1
    bb = dok_array((1, 1))
    bb[0, 0] = 1.0


    # solution from Steve view
    if not args.disable_steve:
        run_solver = not args.dont_solve
        solution = lpsolve(
            A=A,
            B=B,
            b=bb,
            C=C,
            c=cc,
            file_export=args.file_export,
            run_solver=run_solver,
        )
        if run_solver:
            s1, ev1 = solution
            print(s1.tolist())
            print(f"expected value {ev1}")

    # solution from candidate view
    if not args.disable_candidate:
        run_solver = not args.dont_solve
        solution = lpsolve(
            A=-1 * A.T,
            C=B,
            c=bb,
            B=C,
            b=cc,
            file_export=args.file_export,
            run_solver=run_solver,
        )
        if run_solver:
            s1, ev1 = solution
            for i, node in enumerate(c_nodes):
                print(s2.T[i], node)
            print(f"expected value {-ev2}")


def lpsolve(
    A,
    B,
    b,
    C,
    c,
    file_export=None,
    run_solver=True,
):
    # minmax x'Ay st Bx = b, Cy = c and x,y >= 0
    # translate to
    # max c'z st Bx = b, C'z <= A'x

    assert B.shape[0] == b.shape[0]
    assert b.shape[1] == 1

    assert C.shape[1] == A.shape[1]
    assert C.shape[0] == c.shape[0]
    assert c.shape[1] == 1

    model = mathopt.Model()

    # vars
    xvars = [model.add_variable(name=f"x_{i}", lb=0) for i in range(A.shape[0])]
    zvars = [model.add_variable(name=f"z_{i}") for i in range(c.shape[0])]

    # Bx = b
    cons_b = defaultdict(list)
    for (i, j), v in B.items():
        cons_b[i].append(v * xvars[j])
    for i, v in cons_b.items():
        model.add_linear_constraint(sum(v) == b[i, 0])

    # C'z <= A'x
    dA = defaultdict(list)
    dC = defaultdict(list)

    for (j, i), v in C.items():
        dC[i].append(v * zvars[j])
    for (j, i), v in A.items():
        dA[i].append(v * xvars[j])
    for i in range(C.shape[1]):
        eA = dA[i]
        eC = dC[i]
        assert len(eA) >= 1
        assert len(eC) >= 1
        model.add_linear_constraint(sum(eC) <= sum(eA))

    model.maximize(sum(c[i, 0] * zv for (i, zv) in enumerate(zvars)))

    if file_export:
        model_proto = model.export_model()
        model_mps = mps_converter.model_proto_to_mps(model_proto)
        with open(file_export, "w") as f:
            f.write(model_mps)

    if run_solver:
        solver = mathopt.SolverType.GLOP
        result = mathopt.solve(
            model,
            solver,
        )
        if result.termination.reason not in (
            mathopt.TerminationReason.OPTIMAL,
            # mathopt.TerminationReason.FEASIBLE,
        ):
            raise RuntimeError(f"model failed to solve: {result.termination}")

        vars = result.variable_values()
        start = np.array([[vars[v] for v in xvars]])
        optimisation_goal = result.objective_value()
        ev = optimisation_goal
        return (start, ev)


def parse_args():
    parser = argparse.ArgumentParser(
                    prog='Steve game',
                    description='Solves adverserial bisection search game')
    
    default_n=3
    parser.add_argument('-n', '--count', default=default_n,  type=int, help=f"provide problem size, default = {default_n}")

    parser.add_argument('-f', '--file-export', default=None , help="export problem in provided mps file - combine with flags -disable-steve/--disable-candidate")
    parser.add_argument("-b", "--disable-steve", action="store_true" , help="when enabled, model for Steve strategy is not generated" )
    parser.add_argument("-c", "--disable-candidate", action="store_true" , help="when enabled, model for candidate strategy is not generated" )
    parser.add_argument("-s", "--dont-solve", action="store_true" ,help="when enabled, solver does not run, can be used to just export model in mps file"  )

    args = parser.parse_args()
    # validate args
    if  args.count < 1:
        parser.error("flag --count should be >= 1")
    if args.file_export and (args.disable_steve == args.disable_candidate):
        parser.error("excatly one flag of --disable-steve/--disable-candidate should be enabled when exporting to file")

    solve(args)




if __name__ == "__main__":
    parse_args()
