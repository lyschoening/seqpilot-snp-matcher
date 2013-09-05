import argparse
import csv
import subprocess
import numpy
import os
import re
from jinja2 import Environment

__author__ = 'lyschoening'


def read_reference_snps(filename):
    with open(filename, 'r') as snpfile:
        for row in csv.reader(snpfile, delimiter='\t'):
            accession, reference = row[0], row[1]
            if len(reference) > 2:
                reference = '-'
            yield accession, reference

pattern = re.compile(r'.*->\s+([TAGC]+)\s+\((het|homo?)\).*\s(rs\d+)\s.*')



def table_get_snps(file):
    snp_dct = {}
    with open(file, 'r') as table:
        for line in table:
            matches = re.match(pattern, line)

            if matches:
                alt, het_hom, rs_nr = matches.groups()
                snp_dct[rs_nr] = (alt, het_hom)
    return snp_dct





LATEX_SUBS = (
    (re.compile(r'\\'), r'\\textbackslash'),
    (re.compile(r'([{}_#%&$])'), r'\\\1'),
    (re.compile(r'~'), r'\~{}'),
    (re.compile(r'\^'), r'\^{}'),
    (re.compile(r'"'), r"''"),
    (re.compile(r'\.\.\.+'), r'\\ldots'),
)

def escape_tex(value):
    newval = value
    for pattern, replacement in LATEX_SUBS:
        newval = pattern.sub(replacement, newval)
    return newval

def get_template():
    texenv = Environment()
    #    texenv.block_start_string = '((*'
    #    texenv.block_end_string = '*))'
    texenv.variable_start_string = '((('
    texenv.variable_end_string = ')))'
    #    texenv.comment_start_string = '((='
    #    texenv.comment_end_string = '=))'
    texenv.filters['escape_tex'] = escape_tex

    template = texenv.from_string(r"""
\documentclass[11pt]{article}
\usepackage{geometry}                % See geometry.pdf to learn the layout options. There are lots.
\geometry{letterpaper}                   % ... or a4paper or a5paper or ...
\geometry{landscape}                % Activate for for rotated page geometry
%\usepackage[parfill]{parskip}    % Activate to begin paragraphs with an empty line rather than an indent

\usepackage{graphicx}
\usepackage{amssymb}
\usepackage{epstopdf}
\usepackage{longtable, booktabs, multirow}
\heavyrulewidth=.13em
\lightrulewidth=.08em
\cmidrulewidth=.03em

\usepackage{xcolor}
\usepackage{colortbl}
\usepackage[small,bf,singlelinecheck=off]{caption}

\begin{document}
\setlength\LTleft{0pt}
\setlength\LTright{0pt}
\setlength\tabcolsep{2pt}


\newpage

{
\sf
\begin{longtable}{@{\extracolsep{\fill}} l{% for snp in snps %}c{% endfor %}@{\extracolsep{\fill}}}
\toprule
&
{% for accession, reference in snps %}
\rotatebox{85}{((( accession|escape_tex )))} {% if not loop.last %}&{% endif %}
{% endfor %} \\

\cmidrule{2-((( snps|length + 1 ))) }

ID &
{% for accession, reference in snps %}
((( reference|escape_tex ))){% if not loop.last %}&{% endif %}
{% endfor %} \\
\toprule


{% for group in grouped_samples %}
{% for id, snps in group %}
((( id|escape_tex ))) &
{% for snp, type in snps %}
{% if type == 'ref' %}
{% elif type == 'het' %}
\cellcolor{gray!32}
{% elif type == 'homo' %}
\cellcolor{gray!128}
\color{white}
{% endif %}
((( snp )))
{% if not loop.last %}&{% endif %}
{% endfor %} \\

{% endfor %}
{% if not loop.last %}\midrule{% endif %}
{% endfor %}

\bottomrule

\caption{Reference mismatches. Cells in light gray refer to heterozygous variants, cells in dark grey refer to homozygous variants.}
\end{longtable}
}



\newpage

{
\sf
\begin{longtable}{r{% for id, snps in samples %}r{% endfor %}}
\toprule
&


{% for id, snps in samples %}
\rotatebox{90}{\textbf{((( id|escape_tex )))}}
{% if not loop.last %}&{% endif %}
{% endfor %}\\
\toprule

{% for id, snps in samples %}
\textbf{((( id|escape_tex )))} &

{% set i = loop.index0 %}
{% for id, snps in samples %}
{% set j = loop.index0 %}

{% if i == j %}
\multicolumn{1}{c}{--}
{% else %}

{% if mismatches[i][j] == 0 %}
\cellcolor{green!64}
{% else %}
\cellcolor{gray!((( (mismatches[i][j] / snps|length * 50 )|int  )))}
{% endif %}

((( mismatches[i][j]|int ))){% endif %}

{% if not loop.last %}&{% endif %}
{%endfor %} \\
{%endfor %}
\bottomrule

\caption{Sample/sample mismatches (1 per strand).}
\end{longtable}
}

\end{document}
""")

    return template


def main():
    parser = argparse.ArgumentParser(description='Compare a number of JSI Seq Pilot tables')
    parser.add_argument('reference', type=str, help="Reference SNPs, one per line: rs123456<tab>A")
    parser.add_argument('samples', metavar='S', type=str, nargs='+', help="One or more samples to compare")
    parser.add_argument('output', metavar='O', type=str, help="Output prefix")

    args = parser.parse_args()

    snps = list(read_reference_snps(args.reference))
    snps_dct = dict(snps)


    samples = []

    num_samples = len(args.samples)

    mismatches = numpy.zeros((num_samples, num_samples))

    for target in args.samples:
        target_snps = table_get_snps(target)

        def get_sample_list(sample, reference):

            for snp, ref in reference:
                if snp in sample:
                    sample_snp = sample[snp]
                    if sample_snp[1] == 'het':
                        yield (sample_snp[0], 'het')
                    else:
                        yield (sample_snp[0], 'homo')
                else:
                    yield (ref, 'ref')


        samples.append((os.path.splitext(os.path.basename(target))[0], list(get_sample_list(target_snps, snps))))

    def count_mismatches(a, b):
        mismatches = 0
        for a_snp, b_snp in zip(a[1], b[1]):
            if a_snp[0] != b_snp[0]:
                mismatches += 2
            elif a_snp[1] != b_snp[1]:
                mismatches += 1
        return mismatches

    for i, a_sample in enumerate(samples):
        for j, b_sample in enumerate(samples):
            if i == j: continue

            print count_mismatches(a_sample, b_sample)
            mismatches[i][j] = count_mismatches(a_sample, b_sample)

    print mismatches

    print samples


    grouped_samples = []
    taken = []

    for i, sample in enumerate(samples):
        if sample in taken:
            continue

        group = [sample]
        group += [sample for j, sample in enumerate(samples) if i != j and mismatches[i][j] == 0 and sample not in taken]

        for groupitem in group:
            taken.append(groupitem)

        grouped_samples.append(group)

    template = get_template()
    tex_file_name = os.path.abspath("{0}.tex".format(args.output))

    with open(tex_file_name, 'w') as tex_file:
        tex_file.write(template.render(snps=snps, samples=samples, grouped_samples=grouped_samples, mismatches=mismatches))

    for i in range(2): # call twice for proper table layout.
        subprocess.call(('pdflatex', '-output-directory=%s' % os.path.dirname(tex_file_name), tex_file_name))


if __name__ == '__main__':
    main()


