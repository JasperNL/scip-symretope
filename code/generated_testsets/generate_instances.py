import sys
import time
import networkx as nx
import pyscipopt as ps
import pathlib
import glob
from generatorutils import *
from multiprocessing import Pool
import threading
from itertools import chain, combinations, product
import math

if __name__ == "__main__":
    funcs = []

    # Flower snark instances.
    pathlib.Path('instances/flowersnark/').mkdir(parents=True, exist_ok=True)
    def generate_flowersnark_edge3colouring(n):
        sys.stdout = None  # No output from subprocesses.
        ts = time.time()
        c = 3
        outputpath = f"instances/flowersnark/flowersnark{n}_{c}.cip"
        if not pathlib.Path(outputpath).exists():
            g = generate_flower_snark(n)
            make_edge_colouring(g, c, outputpath)
        return time.time() - ts

    for n in range(3, 50, 2):
        funcs.append((f"Flower snark edge 3-colouring; {n:3d}", generate_flowersnark_edge3colouring, (n,)))

    # Run everything on multiple cores
    with Pool(12) as p:
        awaits = []
        for name, func, args in funcs:
            awaits.append((name, p.apply_async(func, args)))

        for name, p in awaits:
            p.wait()
            t = p.get()
            print(f"({t:8.3f}s) {name}")
