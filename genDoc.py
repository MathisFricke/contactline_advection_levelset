import subprocess
import shutil
subprocess.call(["doxygen"])

indexFile = open("html/index.html", "r")
outFile = open("temp.html", "w")

for line in indexFile:
    if (line == "<p>FOR_VIDEO_SEE_README.HTML</p>\n"):
        outFile.write("<center><video src=\"images/animation.ogv\" controls></video></center>")
    else:
        outFile.write(line)
        
indexFile.close()
outFile.close()
shutil.move("temp.html", "html/index.html")

