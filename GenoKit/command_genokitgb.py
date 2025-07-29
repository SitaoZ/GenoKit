# -*- coding: utf-8 -*-
import sys
from GenoKit.version import __version__

def sub_usage(args):
    if len(args) == 1:
        print("          \033[1;33;40m\nUsage  :\033[1m\033[1;35;40m  %s\033[1m" %s (args[0]) )
    elif len(args) == 2:
        print("          \033[1;35;40m%s\033[0m    \033[1;32;40m%s\033[0m" % (args[0], args[1]) )

def main_usage():
    print("\n\033[1;33;40mProgram: \033[0m\033[1;35;40m genBankExtract \033[1;31;40m(pipeline for genome feature extract)\033[0m")
    print("\033[1;33;40mVersion: \033[0m\033[1;32;40m %s\033[0m"%(__version__))
    print("\033[1;33;40mContact: \033[0m\033[1;32;40m Sitao Zhu <zhusitao1990@163.com>\033[0m")
    print("\033[1;33;40mUsage  : \033[0m\033[1;35;40m genBankExtract\033[0m \033[1;31;40m<command> [parameters]\033[0m")
    print("\033[1;33;40mCommand: \033[0m")
    sub_usage(["promoter", "extract promoter from GenBank"])
    sub_usage(["gene    ", "extract gene sequence from GenBank"])
    sub_usage(["rRNA    ", "extract rRNA sequence from GenBank"])
    sub_usage(["tRNA    ", "extract tRNA sequence from GenBank"])
    sub_usage(["UTR     ", "extract UTR  from GenBank"])
    sub_usage(["uORF    ", "extract uORF from GenBank"])
    sub_usage(["CDS     ", "extract CDS  from GenBank"])
    sub_usage(["dORF    ", "extract dORF from GenBank"])
    sub_usage(["exon    ", "extract exon from GenBank"])
    sub_usage(["intron  ", "extract intron from GenBank"])

    sys.exit(1)

def main():
    if len(sys.argv) == 1:
        main_usage()
    elif len(sys.argv) >= 2:
        if sys.argv[1] in ['create','gene','rRNA','tRNA','promoter','UTR','uORF','CDS','dORF','exon','intron']:
            # import 就执行genokit()
            # 安装后，系统存在GenoKit包
            from GenoKit import genokitgb
            #import genokit
        else:
            main_usage()

if __name__ == "__main__":
    main()
