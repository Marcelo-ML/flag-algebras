import sys
from itertools import product

p_modes = [["0.0", "1.0"], ["0.0", "0.5", "1.0"]]

def read_graph_from_file(filename):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]  # remove blank lines

    # Expect "(" as first token and ")" as last
    if lines[0][0] != "(" or lines[-1] != ")":
        raise ValueError("Input must start with '(' and end with ')'.")

    # Read number of vertices
    n = int(lines[0][1])

    edges = []
    # Lines between index 1 and -1 contain edges
    for line in lines[1:-1]:
        a, b = map(int, line.split())
        if not (0 <= a < n and 0 <= b < n):
            raise ValueError(f"Edge ({a}, {b}) contains invalid vertex indices.")
        if a == b:
            raise ValueError(f"Self-loop detected: ({a}, {b})")
        edges.append((a, b))

    return n, edges

def get_useless_masks(n, edges):
    size = 2 ** n
    useless = []

    for a, b in edges:
        # Iterate over all subsets containing both a and b
        # A subset contains a and b if the bits a and b are set: subset & (1<<a) and subset & (1<<b)
        for subset in range(size):
            if (subset & (1 << a)) and (subset & (1 << b)):
                useless.append(subset)

    return useless


def main():
    if len(sys.argv) != 4:
        print("Usage: python gen_graph.py <graph file> <mode> <output file>")
        sys.exit(1)

    graph_path = sys.argv[1]
    n, edges = read_graph_from_file(graph_path)

    graph = ""
    with open(graph_path, "r") as f:
        graph = f.read()

    mode = int(sys.argv[2])
    p_vals = p_modes[mode]

    ## There are values of p that are useless
    ## Let's filter the values that do not yield triangles immediately
    useless_masks = get_useless_masks(n,edges)
    print(useless_masks)

    out_path = sys.argv[3]

    with open(out_path, "a") as file:
        for values in product(p_vals, repeat=(2**n)):
            is_ok = True
            for mask in useless_masks:
                if (values[mask] != "0.0"):
                    is_ok = False
            if (not is_ok):
                continue
            file.write(graph)
            file.write(" ".join(map(str, values)) + "\n")

if __name__ == "__main__":
    main()






















