import argparse

from cargo_oligo_creator import guide, oligo_creator, split_guides_finder

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs="+", type=str, help="Sequences of gRNAs")
    parser.add_argument("-r", "--raw", action="store_true")
    return parser.parse_args()

def main(args):
    split_guides = split_guides_finder.SplitGuidesFinder([guide.Guide(x) for x in args.input]).find_split_guides()
    creator = oligo_creator.OligoCreator(split_guides)
    displays = creator.create_oligos()
    for display in displays:
        print(display.name())
        if args.raw:
            print(display.forward_oligo.display_raw())
            print(display.reverse_oligo.display_raw())
        else:
            print(display.forward_oligo.display_formatted())
            print(display.reverse_oligo.display_formatted())

args = get_args()

if __name__ == "__main__":
    main(args)
