import argparse

from cargo_oligo_creator import guide, oligo_creator, split_guides_finder

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs="+", type=str, help="Sequences of gRNAs")
    return parser.parse_args()

def main(args):
    split_guides = split_guides_finder.SplitGuidesFinder([guide.Guide(x) for x in args.input]).find_split_guides()
    creator = oligo_creator.OligoCreator(split_guides)
    displays = creator.create_oligos()
    for display in displays:
        print(display.name())
        print(display.forward_oligo)
        print(display.reverse_oligo)

args = get_args()

if __name__ == "__main__":
    main(args)
