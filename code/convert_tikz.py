import os
import subprocess


def rem_files(filename, ext_list):
    for i in ext_list:
        try:
            os.remove(filename + i)
        except FileNotFoundError:
            pass


def compile_files(s, filename, rem_file_list, remove_files):
    with open(filename + '.tex', 'w', encoding='utf-8') as f:
        f.write(s)              
    subprocess.call(['latexmk', '-pdf', '-enable-write18', filename + '.tex'])  
    if remove_files:
        rem_files(filename, rem_file_list)


def latex_wrapper(filename):
    texp1 = '\\documentclass{standalone}\n'
    texp2 = '''\\usepackage{tikz}\n\\usepackage{pgfplots}\n\\usepackage{gensymb}
    \n\\usepackage{numprint}\n\\usepackage{amsmath}\n
    \\usetikzlibrary{decorations.pathreplacing}\n
    \\DeclareMathOperator{\\mse}{MSE}\n
    \\tikzset{
    mybrace/.style={decorate,decoration={brace,aspect=#1}}
    }
    \\pgfplotsset{compat=1.16}
    \\usetikzlibrary{backgrounds,external,shapes}
    \\usetikzlibrary{decorations,arrows,calc,arrows.meta,fit,positioning}\n'''

    texp4 = '''\\definecolor{bgColor2}{HTML}{F5793A}
    \\definecolor{dataColor}{HTML}{F5793A}
    \\definecolor{axisColor}{RGB}{0,0,0}
    \\definecolor{labelColor}{gray}{0.30}
    \\definecolor{tickColor}{gray}{0.30}
    \\definecolor{shadeColor}{gray}{0.90}
    \\definecolor{textColor}{RGB}{0,0,0}
    \\definecolor{frameColor}{RGB}{255,255,255}
    \\definecolor{lineColor}{HTML}{0F2080}
    \\definecolor{fillColor}{RGB}{200,200,200}
    \\definecolor{outlineColor}{gray}{0.30}
    \\definecolor{bgColor}{RGB}{200,200,200}
    \\definecolor{reflineColor}{RGB}{220,220,220}
    \\definecolor{lcolor1}{RGB}{216,27,96}
    \\definecolor{lcolor2}{RGB}{30,136,229}
    \\definecolor{lcolor3}{RGB}{255,193,7}
    \\definecolor{lcolor4}{RGB}{0,77,64}
    \\newcommand{\\np}{\\numprint}
    \\begin{document}'''

    latex_str_gl = texp1 + texp2 + texp4
    latex_str_gl += ('\n\\input{"' + filename + '.tex"}')
    latex_str_gl += '\n\\end{document}'
    return latex_str_gl


def create_figures(filename, output_folder, remove_files):
    cwd = os.getcwd()
    os.chdir(output_folder)
    tikz_ext = ['.dpth', '.md5', '.log']
    tex_ext = ['.aux', '.auxlock', '.fdb_latexmk', '.fls', '.log', '.tex']
    latex_str = latex_wrapper(filename)
    tempfigname = 'tempfigs'
    compile_files(latex_str, tempfigname, tex_ext, remove_files)
    try:
        os.remove(filename + '.pdf')
    except FileNotFoundError:
        pass
    os.rename(tempfigname + '.pdf', filename + '.pdf')
    if remove_files:
        rem_files(filename, tikz_ext)
    os.chdir(cwd)


def main():
    create_figures(args.filename, os.path.join(os.getcwd(), args.dir),
                   args.remove_files)
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Convert tikz')
    parser.add_argument('filename')
    parser.add_argument('--dir', default='figures')
    parser.add_argument('--remove_files', default=True, type=bool)
    args = parser.parse_args()
    main()
