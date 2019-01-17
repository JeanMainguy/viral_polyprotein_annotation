import sys


def convert_color(l):
    term_reset = '\x1b[0m'
    html_end_balise = "</span>"
    for col, int in color_dict.items():
        l = l.replace(f'\x1b[1;{int}m', f'<span style="color:{col};font-weight:bold;">')
    l = l.replace(term_reset, html_end_balise)

    return l


if __name__ == '__main__':
    color_dict = {'red': 31,
                  'green': 32,
                  'yellow': 33,
                  'blue': 34,
                  'purple': 35,
                  'teal': 36,
                  'white': 37}
    aln_file = sys.argv[1]
    html_file = aln_file[:-len(".aln")] + '.html'

    html_header = """<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<title>stdin</title>
</head>
<body>
<pre>
"""
    with open(html_file, 'w') as html_wr:
        html_wr.write(html_header)
        for l in open(aln_file):
            html_wr.write(convert_color(l))
        html_wr.write('</pre>')

#
