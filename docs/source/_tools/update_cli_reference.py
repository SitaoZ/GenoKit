"""Generate committed RST options from the current parser, without command dispatch."""
import argparse
import ast
from pathlib import Path
import sys

def load_parser(repo):
    path = repo / 'GenoKit/genokit.py'
    tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    statements = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'args' for t in node.targets):
            break
        statements.append(node)
    else:
        raise RuntimeError('Parser dispatch boundary not found; review generator for the new source layout')
    namespace = {'__name__': '_genokit_docs_parser', '__file__': str(path)}
    previous = sys.argv
    try:
        sys.argv = ['GenoKit', '_documentation_', '_no_dispatch_']
        exec(compile(ast.Module(body=statements, type_ignores=[]), str(path), 'exec'), namespace)
    finally:
        sys.argv = previous
    return namespace['parser'], namespace['_BARE_COMMAND_HELP_PARSERS']

def main():
    options = argparse.ArgumentParser()
    options.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[3])
    options.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1] / 'command_reference.rst')
    args = options.parse_args()
    _, commands = load_parser(args.repo.resolve())
    lines = ['Command reference', '=================', '',
             'This snapshot is generated from the current CLI parser. It includes all',
             '26 public subcommands, their flags, choices and parser defaults. It does',
             'not import command implementations or execute an analysis. Defaults of',
             '``None`` mean unspecified; some commands require an explicit output at',
             'runtime. See the user-guide pages for biological semantics and known',
             'implementation limits; parser acceptance is not a guarantee of behavior.', '',
             'Regenerate after changing the CLI using :doc:`documentation`.', '']
    for name, parser in commands.items():
        title = 'GenoKit ' + name
        lines += ['.. _cli-' + name + ':', '', title, '-' * len(title), '', '.. code-block:: text', '']
        parser.epilog = None  # Avoid repeating historical epilog examples as normative instructions.
        lines.extend('   ' + line for line in parser.format_help().strip().splitlines())
        lines += ['', '**Parser defaults and requirements**', '']
        for action in parser._actions:
            if isinstance(action, argparse._HelpAction):
                continue
            flags = ', '.join(action.option_strings)
            default = 'required' if action.required else 'default: ' + repr(action.default)
            choices = '; choices: ' + ', '.join(map(str, action.choices)) if action.choices is not None else ''
            lines.append('* ``' + flags + '``: ' + default.replace('None', 'unspecified') + choices + '.')
        lines.append('')
    args.output.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print('Generated', len(commands), 'commands:', args.output)
if __name__ == '__main__':
    main()
