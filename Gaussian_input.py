input_file = "BR_cart.gau"
settings = "Gaussian_Input1.txt"
charge = "0"
multiplicity = "1"
c_m = charge + " " + multiplicity

n = -1
with open("%s" %input_file, "r") as my_file:
    for line in my_file:
        if "BR" in line:
            words = line.split()
            title = words[2]
            n += 1
            with open("%s" % settings, "r") as f:
                with open("Water_R_%s.com" % n, "w") as f2:
                    for line1 in f:
                        f2.write(line1)
                    f2.write("Water with R = %s Angstroms\n\n" % title)
                    f2.write(c_m + "\n")
        elif "H" or "O" in line:
            with open("Water_R_%s.com" % n, "a+") as f2:
                f2.write(line)