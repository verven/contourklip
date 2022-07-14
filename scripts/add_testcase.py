import re
import pathlib
import argparse

def svg_prefix(viewBox):

    return f'<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n<svg width="100%" height="100%" viewBox="{viewBox}" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" style="fill-rule:evenodd;clip-rule:evenodd;stroke-linecap:round;">\n'

def svg_suffix():
    return "</svg>"



def to_testcase_svg(svg_str):
    svg_str = svg_str.strip()
    path_re = "<path(?:.*?)d=\"(.*?)\"(?:.*?)\/>"
    path_re= re.compile(path_re)
    viewbox_re = r"svg(?:.*?)viewBox=\"(.*?)\"(?:.*?)"
    viewbox_re = re.compile(viewbox_re)

    viewbox_vals = re.search(viewbox_re, svg_str)
    if not viewbox_vals:
        print("couldn't infer svg dimensions")
        exit(1)
    t= viewbox_vals.group(1).strip().split()[-2:]
    height, width = int(t[0]), int(t[1])

    out= svg_prefix(viewbox_vals.group(1))

    for idx, p in enumerate(re.finditer(path_re, svg_str)):
        print(p.group(1), "\n\n")
        out += f'<path d="{p.group(1)}" style="fill:lightgrey;opacity:0.5;stroke:black;stroke-width:{min(height, width)/1000}px;"/>\n'
        # if idx ==1:
        #     break

    return out +svg_suffix()



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process an svg file so that it creates a testcase svg file.')
    parser.add_argument('filepath', metavar='file', type=str, help='the filepath to the svg file')
    args = parser.parse_args()
    p= args.filepath

    if not p.endswith(".svg"):
        print("the file is not an .svg file")
        exit(1)
    out_dir = pathlib.Path(p[:-4])
    out_dir.mkdir(parents=True, exist_ok=True)
    write_p = out_dir / "in.svg"
    with open(p) as file_in:
        svg_str = file_in.read()
        with write_p.open("w", encoding="utf-8") as f:
            f.write(to_testcase_svg(svg_str))
