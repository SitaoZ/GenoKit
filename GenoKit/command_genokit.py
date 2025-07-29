#!/usr/bin/env python3
import argparse
import sys

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
          colorize(" v0.2.6", "LIGHT_YELLOW"))
    print(colorize("  Contact: Sitao Zhu <zhusitao1990@163.com>", "LIGHT_CYAN"))
    print()

    # 使用说明
    print(colorize("  Usage:", "LIGHT_MAGENTA", style="BOLD"))
    print("    " + colorize("GenoKit <command> [parameters]", "LIGHT_WHITE"))
    print()

    # 命令分类
    sections = [
        ("Database", [
            ("create", "Create GFF/GTF database"),
            ("stat", "Database statistics")
        ]),
        ("Extract", [
            ("gene", "Extract gene sequence"),
            ("mrna", "Extract mRNA sequence"),
            ("transcript", "Extract transcript sequence"),
            ("exon", "Extract exon sequence"),
            ("intron", "Extract intron sequence"),
            ("cds", "Extract CDS sequence"),
            ("utr", "Extract 5'/3'UTR sequence"),
            ("uorf", "Extract uORF sequence"),
            ("dorf", "Extract dORF sequence"),
            ("promoter", "Extract promoter sequence"),
            ("terminator", "Extract terminator sequence"),
            ("igr", "Extract intergenic region")
        ]),
        ("Design", [
            ("primer", "Primer design"),
            ("sgrna", "Single guide RNA design"),
            ("sirna", "Small interfering RNA design"),
            ("motif", "Motif search")
        ]),
        ("Visualize", [
            ("vision", "Snapshot gene structure"),
            ("circos", "Circlize genome structure")
        ])
    ]

    # 打印命令
    for section, commands in sections:
        print(colorize(f"  {section}:", "LIGHT_GREEN", style="BOLD"))
        for cmd, desc in commands:
            print("    " +
                  colorize(f"{cmd.ljust(12)}", "LIGHT_BLUE", style="BOLD") +
                  colorize(desc, "WHITE"))
        print()

    # 标志
    print(colorize("  Flags:", "LIGHT_MAGENTA", style="BOLD"))
    print("    " +
          colorize("-h, --help".ljust(12), "LIGHT_YELLOW", style="BOLD") +
          colorize("Show this help message", "WHITE"))
    print()

def main():
    if len(sys.argv) == 1:
        main_usage()
        sys.exit(0)
    elif len(sys.argv) >= 2:
        if sys.argv[1] in ['create','stat','gene','mrna','transcript','igr','promoter',
                           'sirna', 'sgrna', 'primer', 'motif', 'vision', 'circos',
                           'terminator','utr','uorf','cds','dorf','exon','intron']:
            # import 就执行genokit()
            # 安装后，系统存在GenoKit包
            from GenoKit import genokit
            #import genokit
        else:
            main_usage()

    # 其他命令处理逻辑...
if __name__ == '__main__':
    main()
