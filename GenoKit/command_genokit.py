#!/usr/bin/env python3

import os
import sys
from GenoKit.version import __version__, __date__

class Colors:
    """ANSI颜色代码"""
    # 基础颜色
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'

    # 亮色（更鲜艳）
    LIGHT_RED = '\033[91m'
    LIGHT_GREEN = '\033[92m'
    LIGHT_YELLOW = '\033[93m'
    LIGHT_BLUE = '\033[94m'
    LIGHT_MAGENTA = '\033[95m'
    LIGHT_CYAN = '\033[96m'
    LIGHT_WHITE = '\033[97m'

    # 背景色
    BG_BLUE = '\033[44m'
    BG_CYAN = '\033[46m'

    # 样式
    BOLD = '\033[1m'
    DIM = '\033[2m'
    ITALIC = '\033[3m'
    UNDERLINE = '\033[4m'
    BLINK = '\033[5m'
    REVERSE = '\033[7m'

    # 重置所有样式
    RESET = '\033[0m'

def colorize(text, color=None, bg_color=None, style=None):
    """应用颜色和样式"""
    if not sys.stdout.isatty() or os.environ.get('NO_COLOR') is not None:
        return text
    color_code = getattr(Colors, color.upper(), '') if color else ''
    bg_code = getattr(Colors, bg_color.upper(), '') if bg_color else ''
    style_code = getattr(Colors, style.upper(), '') if style else ''
    return f"{style_code}{bg_code}{color_code}{text}{Colors.RESET}"

def print_banner_1():
    """打印程序横幅"""
    banner = r"""
  ██████  ███████ ███    ██  ██████  ██   ██ ██ ████████ 
 ██       ██      ████   ██ ██    ██ ██  ██  ██    ██    
 ██   ██  █████   ██ ██  ██ ██    ██ █████   ██    ██    
 ██    ██ ██      ██  ██ ██ ██    ██ ██  ██  ██    ██    
  ██████  ███████ ██   ████  ██████  ██   ██ ██    ██    
"""
    print(colorize(banner, "LIGHT_BLUE", style="BOLD"))

def print_banner_2():
    """打印程序横幅"""
    banner = r"""
  ██████  ███████ ███    ██  ██████  ██   ██ ██ ████████ 
 ██       ██      ████   ██ ██    ██ ██  ██  ██    ██    
 ██   ██  █████   ██ ██  ██ ██    ██ █████   ██    ██    
 ██    ██ ██      ██  ██ ██ ██    ██ ██  ██  ██    ██    
  ██████  ███████ ██   ████  ██████  ██   ██ ██    ██    
    """
    print(colorize(banner, "LIGHT_WHITE", style="BOLD"))

def print_banner():
    banner = r"""
   ______                 __ __ _ __ 
  / ____/__  ____  ____  / //_/(_) /_
 / / __/ _ \/ __ \/ __ \/ ,<  / / __/
/ /_/ /  __/ / / / /_/ / /| |/ / /_  
\____/\___/_/ /_/\____/_/ |_/_/\__/ 
    """
    print(colorize(banner, "LIGHT_BLUE", style="BOLD"))

def print_banner_standard():
    banner = r"""
   ____                  _  ___ _   
  / ___| ___ _ __   ___ | |/ (_) |_ 
 | |  _ / _ \ '_ \ / _ \| ' /| | __|
 | |_| |  __/ | | | (_) | . \| | |_ 
  \____|\___|_| |_|\___/|_|\_\_|\__|
                                    
    """
    print(colorize(banner, "LIGHT_BLUE", style="BOLD"))

def print_banner_ivrit():
    banner = r"""
   ____                  _  ___ _   
  / ___| ___ _ __   ___ | |/ (_) |_ 
 | |  _ / _ \ '_ \ / _ \| ' /| | __|
 | |_| |  __/ | | | (_) | . \| | |_ 
  \____|\___|_| |_|\___/|_|\_\_|\__|
                                    

"""
    print(colorize(banner, "LIGHT_BLUE", style="BOLD"))

def print_banner_shadow():
    banner = r"""
   ___|                      |  / _)  |   
  |       _ \  __ \    _ \   ' /   |  __| 
  |   |   __/  |   |  (   |  . \   |  |   
 \____| \___| _|  _| \___/  _|\_\ _| \__| 
                                          
    """
    print(colorize(banner, "LIGHT_BLUE", style="BOLD"))

def print_banner_stop():
    banner = r"""
  ______                   _    _ _      
 / _____)                 | |  / |_)_    
| /  ___  ____ ____   ___ | | / / _| |_  
| | (___)/ _  )  _ \ / _ \| |< < | |  _) 
| \____/( (/ /| | | | |_| | | \ \| | |__ 
 \_____/ \____)_| |_|\___/|_|  \_)_|\___)
                                         
    """
    print(colorize(banner, "LIGHT_BLUE", style="BOLD"))
    
def print_banner_3d():
    banner = r"""
  ██████╗ ███████╗███╗   ██╗ ██████╗ ██╗  ██╗██╗████████╗
 ██╔════╝ ██╔════╝████╗  ██║██╔═══██╗██║ ██╔╝██║╚══██╔══╝
 ██║  ██╗ █████╗  ██╔██╗ ██║██║   ██║█████╔╝ ██║   ██║   
 ██║  ╚██╗██╔══╝  ██║╚██╗██║██║   ██║██╔═██╗ ██║   ██║   
 ╚██████╔╝███████╗██║ ╚████║╚██████╔╝██║  ██╗██║   ██║   
  ╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝ ╚═╝  ╚═╝╚═╝   ╚═╝   
    """
    print(colorize(banner, "LIGHT_WHITE", style="BOLD"))

def main_usage():
    """打印带颜色的帮助信息"""
    # 打印横幅
    print_banner()
    # print_banner_standard()
    # print_banner_shadow()
    # print_banner_stop()

    # 版本和联系信息
    print(colorize("  GenoKit - Genome and Gene ToolKit", "LIGHT_CYAN", style="BOLD") +
          colorize(f"  v{__version__}", "LIGHT_YELLOW"))
    print(colorize("  Contact: Sitao Zhu <zhusitao1990@163.com>", "LIGHT_CYAN"))
    print()

    # 使用说明
    print(colorize("  Usage:", "LIGHT_MAGENTA", style="BOLD"))
    print("    " + colorize("GenoKit [--help | --version]", "LIGHT_WHITE"))
    print("    " + colorize("GenoKit <command>  [options]", "LIGHT_WHITE"))
    print()

    # 命令分类
    sections = [
        ("Database", [
            ("create", "Build a genome annotation database"),
            ("fasta", "Extract selected FASTA records"),
            ("iupac", "Display standard IUPAC codes"),
            ("kmer", "Build genome/transcriptome k-mer catalogues"),
            ("stat", "Summarize annotation statistics")
        ]),
        ("Extract", [
            ("gene", "Extract genomic gene sequences"),
            # ("mrna", "Extract messenger RNA sequences"),
            ("transcript", "Extract spliced transcript sequences"),
            ("exon", "Extract transcript exon sequences"),
            ("intron", "Extract gene intron sequences"),
            ("cds", "Extract coding sequence regions"),
            ("utr", "Extract 5′ and 3′ UTR sequences"),
            ("uorf", "Extract upstream open reading frames"),
            ("dorf", "Extract downstream open reading frames"),
            ("promoter", "Extract promoter sequence regions"),
            ("terminator", "Extract terminator sequence regions"),
            ("igr", "Extract intergenic region sequences")
        ]),
        ("Design", [
            ("barcode", "Design DNA and RNA barcode sets"),
            ("motif", "Search and analyze sequence motifs"),
            ("primer", "Design gene-specific primer pairs"),
            ("sgrna", "Design CRISPR single-guide RNAs"),
            ("shrna", "Design small hairpin RNA candidates"),
            ("sirna", "Design small interfering RNAs")
        ]),
        ("Visualize", [
            ("circos", "Visualize circular genome structure"),
            ("hgvs", "Visualize DNA, RNA, and protein variants"),
            ("ppi", "Visualize protein interaction networks"),
            ("view", "Visualize annotated gene structure"),
            # ("expression", "Gene expression pattern")
        ])
    ]

    # 打印命令
    for section, commands in sections:
        print(colorize(f"  {section}:", "LIGHT_GREEN", style="BOLD"))
        for cmd, desc in commands:
            print("    " +
                  colorize(f"{cmd.ljust(15)}", "LIGHT_BLUE", style="BOLD") +
                  colorize(desc, "WHITE"))
        # print()

    # 标志
    print()
    print(colorize("  Flags:", "LIGHT_MAGENTA", style="BOLD"))
    print("    " +
          colorize("-h, --help".ljust(15), "LIGHT_YELLOW", style="BOLD") +
          colorize("Display command help information", "WHITE"))
    print("    " +
          colorize("-v, --version".ljust(15), "LIGHT_YELLOW", style="BOLD") +
          colorize("Display GenoKit version information", "WHITE"))
    print()
    print("  Run 'GenoKit <command> -h' for details on a command.")
    # print(colorize("  Examples:", "LIGHT_MAGENTA", style="BOLD"))
    # print("    GenoKit create -h")
    # print("    GenoKit create -s gtf -a human.gtf -g human.fa -o ./ -p human")
    # print()


def main():

    COMMANDS = ['create', 'fasta', 'iupac', 'kmer', 'stat',
                'gene', 'transcript', 'igr', 'promoter', 'terminator',
                'utr', 'uorf', 'cds', 'dorf', 'exon', 'intron',
                'barcode', 'motif', 'primer', 'sgrna', 'shrna', 'sirna',
                'circos', 'hgvs', 'ppi', 'view'
                ]
    if len(sys.argv) == 1:
        main_usage()
        return

    cmd = sys.argv[1]
    if cmd in COMMANDS:
        # Importing genokit parses and dispatches the selected subcommand.
        from GenoKit import genokit
        return
    if cmd in ('-h', '--help'):
        main_usage()
        return
    if cmd in ('-v', '--version'):
        print("GenoKit v{} ({})".format(__version__, __date__))
        return

    print(f"error: unknown command '{cmd}'", file=sys.stderr)
    print("Run 'GenoKit -h' to see available extract.", file=sys.stderr)
    raise SystemExit(2)
if __name__ == '__main__':
    main()
