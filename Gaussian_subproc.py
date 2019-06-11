import subprocess as proc
for i in range(13):
    proc.call(["g09", "Water_R_%s.com" %i])