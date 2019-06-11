with open("SPE_W_data.txt", "w") as my_file:
    for i in range(13):
        with open("Water_R_%s.log" %i, "r") as f:
            for line in f:
                words = line.split()
                if len(words) > 0 and words[0] == "Water":
                    BR = words[4]
                    my_file.write(BR + " ")
                elif "E(RHF)" in line:
                    SPE = words[4]
                    my_file.write(SPE + "\n")