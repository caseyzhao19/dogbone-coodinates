import os

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# DO NOT RUN UNTIL YOU HAVE MADE A COPY OF PLOTS SO THAT IF THE RENAMING IS BOTCHED YOU DON'T LOSE ALL YOUR WORK
# ALSO, REMEMBER TO COPY THE PRINTOUT FROM DOG BONE TO YOUR TEXT DOCUMENT
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

for file in os.listdir("./plots"):
    num = file[4:-4]
    if 0 <= int(num) <= 9:
        os.rename("./plots/" + file, "./plots/plot00" + num + ".jpg")
    elif 10 <= int(num) <= 99:
        os.rename("./plots/" + file, "./plots/plot0" + num + ".jpg")
