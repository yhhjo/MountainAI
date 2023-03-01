import numpy as np
import argparse 
import SampleGenerator
import Orogens

def main():
    parser = argparse.ArgumentParser(description='Generate geothermal and tectonic histories.')
    parser.add_argument('-n', type=int, help='number of samples')
    parser.add_argument('-c', type=int, help='chunk size')
    parser.add_argument('-p', type=str, help=f'orogen {Orogens.OROGENS}')
    parser.add_argument('-o', type=str, help='output directory')
    parser.add_argument('-t', type=str, help='temporary directory')

    args = parser.parse_args()

    if not all([args.n, args.c, args.p, args.o, args.t]):
        parser.print_help()
        exit()

    if args.n and args.n < 1:
        parser.error('Invalid value for n. Must be >1')
    if args.c and args.c < args.n:
        parser.error('Invalid value for c. c < n')
    if args.p not in Orogens.OROGENS:
        parser.error(f'Invalid value for p. Options: {Orogens.OROGENS}')

    n = args.n if args.n else 100
    chunk_size = args.c if args.c else 100
    orogen = args.p if args.p else "CC"
    output_dir = args.o if args.o else "out"
    tmp_dir = args.t if args.t else "tmp"

    SampleGenerator.generate_models(n, chunk_size, orogen, output_dir, tmp_dir)


if __name__ == '__main__':
    main()
