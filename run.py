import argparse

from cargo_oligo_creator import guide, oligo_creator, split_guides_finder, genbank_writer

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs="+", type=str, help="Sequences of gRNAs")
    parser.add_argument("-r", "--raw", action="store_true")
    parser.add_argument("-g", "--genbank", type=str, help="Path to store generated Genbank files")
    return parser.parse_args()

def main(args):
    split_guides = split_guides_finder.SplitGuidesFinder([guide.Guide(x) for x in args.input]).find_split_guides()
    creator = oligo_creator.OligoCreator(split_guides)
    displays = creator.create_oligos()
    if args.genbank is not None:
        writer = genbank_writer.GenbankWriter(args.genbank, displays)
        writer.write()
        if writer.success:
            print("wrote {} files to {}".format(len(displays)*2, args.genbank))
    else:
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
