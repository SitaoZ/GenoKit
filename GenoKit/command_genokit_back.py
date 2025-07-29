# -*- coding: utf-8 -*-
import sys
from GenoKit.version import __version__

def sub_usage(args):
    if len(args) == 1:
        print("          \033[1;33;40m\nUsage  :\033[1m\033[1;35;40m  %s\033[1m" % (args[0]) )
    elif len(args) == 2:
        print("          \033[1;35;40m%s\033[0m    \033[1;32;40m%s\033[0m" % (args[0], args[1]) )

def main_usage():
    print("\n\033[1;33;40mProgram: \033[0m\033[1;35;40m GenoKit \033[1;31;40m(Genome and Gene ToolKit)\033[0m")
    print("\033[1;33;40mVersion: \033[0m\033[1;32;40m %s\033[0m"%(__version__))
    print("\033[1;33;40mContact: \033[0m\033[1;32;40m Sitao Zhu <zhusitao1990@163.com>\033[0m")
    print("\033[1;33;40mUsage  : \033[0m\033[1;35;40m GenoKit\033[0m \033[1;31;40m<command> [parameters] \033[0m")
    print("\033[1;33;40mCommand: \033[0m")
    sub_usage(["create    ", "create GFF/GTF database"])
    sub_usage(["stat      ", "database statistics"])
    sub_usage(["cds       ", "extract CDS sequence"])
    sub_usage(["dorf      ", "extract dORF sequence"])
    sub_usage(["exon      ", "extract exon sequence"])
    sub_usage(["gene      ", "extract gene sequence"])
    sub_usage(["intron    ", "extract intron sequence"])
    sub_usage(["igr       ", "extract intergenic region"])
    sub_usage(["mrna      ", "extract mRNA sequence"])
    sub_usage(["motif     ", "motif search"])
    sub_usage(["primer    ", "primer design"])
    sub_usage(["promoter  ", "extract promoter sequence"])
    sub_usage(["sgrna     ", "single guide rna design"])
    sub_usage(["sirna     ", "small interfering rna design"])
    sub_usage(["terminator", "extract terminator sequence"])
    sub_usage(["transcript", "extract transcript sequence"])
    sub_usage(["uorf      ", "extract uORF sequence"])
    sub_usage(["utr       ", "extract 5/3UTR sequence"])
    sub_usage(["vision    ", "visualize gene structure"])
    # command_genokit.py
    print("\033[1;33;40mFlag: \033[0m")
    sub_usage(["-h, --help", "help for GenoKit"])
    
    sys.exit(1)


def main():
    if len(sys.argv) == 1:
        main_usage()
    elif len(sys.argv) >= 2:
        if sys.argv[1] in ['create','stat','gene','mrna','transcript','igr','promoter',
                           'sirna', 'sgrna', 'primer', 'motif', 'vision',
                           'terminator','utr','uorf','cds','dorf','exon','intron']:
            # import 就执行genokit()
            # 安装后，系统存在GenoKit包
            from GenoKit import genokit
            #import genokit
        else:
            main_usage()

if __name__ == "__main__":
    main()
