import networkx as nx
import pyscipopt as ps
import pathlib
from itertools import product, combinations, permutations


def get_graph(path, output=False):
    """Returns a graph in .col-format as NetworkX graph"""
    nnodes = -1
    nedges = -1
    g : nx.Graph = nx.Graph()

    with open(path, "r") as f:
        for line in f:
            if line.strip() == "":
                continue
            if line[0] == "c":
                # Comment line
                continue
            if line[0] == "p":
                linesl = line.strip().split()
                assert(linesl[1] in ["edge", "edges"])
                nnodes = int(linesl[2])
                nedges = int(linesl[3])
                continue
            if line[0] == "e":
                linesl = line.strip().split()
                i = int(linesl[1])
                j = int(linesl[2])
                g.add_edge(i, j)
                continue
            linesl = line.rstrip("\n").split()
            if output:
                print(linesl)
            raise Exception("Unknown line type", line.strip())

    if output:
        if len(g.nodes) != nnodes:
            print(f"Different number of nodes: {len(g.nodes)} gotten, {nnodes} expected ({pathlib.Path(path).name})")
        if len(g.edges) not in [nedges, nedges/2]:
            print(f"Different number of edges: {len(g.edges)} gotten, {nedges} expected ({pathlib.Path(path).name})")
    # assert(len(g.nodes) == nnodes)
    # assert(len(g.edges) == nedges)
    return g


def generate_flower_snark(n, return_nodesets=False):
    nodes = list(range(4 * n))
    A = nodes[0 * n : 1 * n]
    B = nodes[1 * n : 2 * n]
    C = nodes[2 * n : 3 * n]
    D = nodes[3 * n : 4 * n]

    g : nx.Graph = nx.Graph()
    g.add_nodes_from(nodes)

    # Build n copies of the star graph on 4 vertices, with center A, and surrounding edges B, C, D.
    for i in range(n):
        g.add_edges_from([
            (A[i], B[i]),
            (A[i], C[i]),
            (A[i], D[i])
        ])

    # Construct the n-cycle on the nodes (B1,...,Bn).
    g.add_edges_from([
        (B[i - 1], B[i])
        for i in range(n)
    ])

    # Construct the 2n cycle (C1,...,Cn,D1,...,Dn).
    for i in range(n):
        if i == 0:
            g.add_edge(D[-1], C[0])
            g.add_edge(C[-1], D[0])
        else:
            g.add_edge(C[i - 1], C[i])
            g.add_edge(D[i - 1], D[i])

    if return_nodesets:
        return g, (A, B, C, D)
    else:
        return g


# Utility stuff
def gcd(x, y):

    if x > y:
        small = y
    else:
        small = x
    for i in range(1, small+1):
        if((x % i == 0) and (y % i == 0)):
            gcd = i

    return gcd


def get_generator_cycle_notation(gen):
    s = ""
    seen = set()
    for i in range(len(gen)):
        if i in seen:
            continue
        s += "("
        j = i
        s += str(i)
        j = gen[i]
        while j != i:
            s += ", "
            s += str(j)
            j = gen[j]
        s += ")"
    return s


def make_stable_set(g, outputpath, maxncons=None, maxnvars=None):
    model = ps.Model()

    if maxnvars is not None and len(g.nodes) > maxnvars:
        return
    if maxncons is not None and len(g.edges) > maxncons:
        return

    x = {}
    for i in sorted(g.nodes):
        args = str(i).strip(" ()").replace(" ", "")
        x[i] = model.addVar(name=f"x[{args}]", vtype="B")

    for (i, j) in sorted(g.edges):
        model.addCons(x[i] + x[j] <= 1)

    model.setObjective(ps.quicksum(x.values()), sense="maximize")

    model.writeProblem(outputpath)

    model.freeProb()
    del model


def make_max_kcolourable_subgraph(g, k, outputpath=None, maxncons=None, maxnvars=None, return_model=False, order=None):
    model = ps.Model()

    if maxnvars is not None and len(g.nodes) * k > maxnvars:
        return
    if maxncons is not None and len(g.edges) * k > maxncons:
        return

    if order is None:
        order = sorted(g.nodes)
    else:
        assert all(node in g.nodes for node in order)
        assert len(order) == len(g.nodes)

    colours = range(k)

    x = {}
    # Add the variables for each variable in the suggested order.
    for i in order:
        for c in colours:
            args = str(i).strip(" ()").replace(" ", "")
            args = f"{args},{c}"
            x[i, c] = model.addVar(name=f"x[{args}]", vtype="B")

    for (i, j) in sorted(g.edges):
        for c in colours:
            model.addCons(x[i, c] + x[j, c] <= 1)

    for i in order:
        model.addCons(ps.quicksum(x[i, c] for c in colours) <= 1)

    model.setObjective(ps.quicksum(x.values()), sense="maximize")

    if outputpath is not None:
        model.writeProblem(outputpath)

    if return_model:
        return model, (colours, x)
    else:
        model.freeProb()
        del model


def make_vertex_colouring(g, k, outputpath):
    model = ps.Model()

    colours = range(k)

    x = {}
    for i in sorted(g.nodes):
        for c in colours:
            args = str(i).strip(" ()").replace(" ", "")
            args = f"{args},{c}"
            x[i, c] = model.addVar(name=f"x[{args}]", vtype="B")

    for (i, j) in sorted(g.edges):
        for c in colours:
            model.addCons(x[i, c] + x[j, c] <= 1)

    for i in sorted(g.nodes):
        model.addCons(ps.quicksum(x[i, c] for c in colours) == 1)

    model.setObjective(ps.quicksum(x.values()), sense="maximize")

    model.writeProblem(outputpath)

    # model.freeProb()
    # del model
    return model


def make_edge_colouring(g, c, outputpath=None, maxncons=None, maxnvars=None, return_model=False, nodeorder=None, edgeorder=None):
    model = ps.Model()

    if maxnvars is not None and len(g.nodes) * c > maxnvars:
        return
    if maxncons is not None and len(g.edges) * c > maxncons:
        return

    if nodeorder is None:
        nodeorder = sorted(g.nodes)
    else:
        assert all(node in g.nodes for node in nodeorder)
        assert len(nodeorder) == len(g.nodes)

    if edgeorder is None:
        edgeorder = sorted(g.edges)
    else:
        assert all(edge in g.edges for edge in edgeorder)
        assert len(edgeorder) == len(g.edges)

    colours = range(c)

    x = {}
    for i, j in edgeorder:
        argsi = str(i).strip(" ()").replace(" ", "")
        argsj = str(j).strip(" ()").replace(" ", "")
        for c in colours:
            x[tuple(sorted((i, j))), c] = model.addVar(name=f"x[{argsi},{argsj},{c}]", vtype="B")

    for i in nodeorder:
        for c in colours:
            model.addCons(ps.quicksum(x[tuple(sorted((i, j))), c] for j in g.neighbors(i)) <= 1 )

    for i, j in edgeorder:
        model.addCons(ps.quicksum(x[tuple(sorted((i, j))), c] for c in colours) == 1 )

    # model.setObjective(ps.quicksum(x.values()), sense="maximize")

    if outputpath is not None:
        model.writeProblem(outputpath)

    if return_model:
        return model, (colours, x)
    else:
        model.freeProb()
        del model
