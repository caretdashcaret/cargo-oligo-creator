import argparse

from cargo_oligo_creator import guide, oligo_creator

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs="+", type=str, help="Sequences of gRNAs")
    return parser.parse_args()

def main(args):
    creator = oligo_creator.OligoCreator([guide.Guide(x) for x in args.input])
    oligos = creator.create_oligos()
    for oligo in oligos:
        print(oligo)

args = get_args()

if __name__ == "__main__":
    main(args)
