import glob
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--sector", type=int, required=True)
args = parser.parse_args()

cwd = os.getcwd()
val = f"sector{args.sector}"

flagfiles = [os.path.basename(f) for f in glob.glob("flag/*jpg")]
singlefiles = [os.path.basename(f) for f in glob.glob("single/*jpg")]
multfiles = [os.path.basename(f) for f in glob.glob("multiple/*jpg")]

numflag = len(flagfiles)
numsingle = len(singlefiles)
nummult = len(multfiles)

htmlsummary = f"""<html>
<! -- -->

<head>
<title></title>
<body bgcolor="#ffffff">
<h2>{val} Candidates</h1></body>

<a href={val}multipleplots.html>{nummult} multi-point candidates</a>
<br>
<a href={val}singleplots.html>{numsingle} single-point candidates (less likely)</a>
<br>
<a href={val}flagplots.html>{numflag} flagged candidates</a>

</html>
"""

with open(f"{val}resultssummary.html", "w") as f:
    f.write(htmlsummary)

htmlmult = f"""<html>
<! -- -->

<head>
<title>Multi-point Candidates</title>
<body bgcolor="#ffffff">
<h2>{val}</h1></body>"""
for mf in multfiles:
    htmlmult += f"""<figure>\n"""
    htmlmult += f"""  <img src=multiple/{mf}>\n"""
    htmlmult += f"""  <figcaption>{mf}</figcaption>\n"""
    htmlmult += f"""</figure>\n"""
htmlmult += """</html>"""

htmlsingle = f"""<html>
<! -- -->

<head>
<title>Single Point Candidates</title>
<body bgcolor="#ffffff">
<h2>{val} Single Point</h1></body>"""
for sf in singlefiles:
    htmlsingle += f"""<figure>\n"""
    htmlsingle += f"""  <img src=single/{sf}>\n"""
    htmlsingle += f"""  <figcaption>{sf}</figcaption>\n"""
    htmlsingle += f"""</figure>\n"""
htmlsingle += """</html>"""

htmlflag = f"""<html>
<! -- -->

<head>
<title>Flagged Candidates</title>
<body bgcolor="#ffffff">
<h2>{val} Flagged Candidates</h1></body>"""
for ff in flagfiles:
    htmlflag += f"""<figure>\n"""
    htmlflag += f"""  <img src=flag/{ff}>\n"""
    htmlflag += f"""  <figcaption>{ff}</figcaption>\n"""
    htmlflag += f"""</figure>\n"""
htmlflag += """</html>"""

with open(f"{val}multipleplots.html", "w") as f:
    f.write(htmlmult)
with open(f"{val}singleplots.html", "w") as f:
    f.write(htmlsingle)
with open(f"{val}flagplots.html", "w") as f:
    f.write(htmlflag)
