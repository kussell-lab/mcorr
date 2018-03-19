from __future__ import print_function
from argparse import ArgumentParser

def main():
    """Run fitting using lmfit"""
    parser = ArgumentParser(description="Convert to distant matrix")
    parser.add_argument("fit_res_file", type=str)
    parser.add_argument("output_file", type=str)
    parser.add_argument('--by', nargs='?', const="theta", type=str, default="theta")
    opts = parser.parse_args()
    datafile = opts.fit_res_file
    outfile = opts.output_file
    byvalue = opts.by

    dmap = {}
    with open(datafile) as reader:
        header = reader.readline().rstrip().split(",")
        for line in reader:
            terms = line.rstrip().split(",")
            group = terms[0]
            if "_vs_" in group:
                isolates = group.split("_vs_")
                ddmap = dmap.get(isolates[0], {})
                ddmap[isolates[1]] = terms[header.index(byvalue)]
                dmap[isolates[0]] = ddmap

                ddmap = dmap.get(isolates[1], {})
                ddmap[isolates[0]] = terms[header.index(byvalue)]
                dmap[isolates[1]] = ddmap
    isolates = sorted(dmap.keys())
    with open(outfile, 'w') as writer:
        writer.write("," + ",".join(isolates) + "\n")
        for isolate1 in isolates:
            writer.write(isolate1)
            for isolate2 in isolates:
                if isolate1 == isolate2:
                    value = 0
                else:
                    value = float(dmap[isolate1][isolate2])
                writer.write(",%g" % value)
            writer.write("\n")
            


if __name__ == "__main__":
    main()