import os

for file in os.listdir("./plots"):
    num = file[4:-4]
    if 0 <= int(num) <= 9:
        os.rename("./plots/" + file, "./plots/plot00" + num + ".jpg")
    elif 10 <= int(num) <= 99:
        os.rename("./plots/" + file, "./plots/plot0" + num + ".jpg")
