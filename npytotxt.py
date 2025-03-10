import numpy as np
import os

dirPath = "g64"

n = 1
while(True):
    filedx = f"{dirPath}/{n:05d}_dx.npy"
    filedy = f"{dirPath}/{n:05d}_dy.npy"
    dArrayx = np.load(filedx)
    dArrayy = np.load(filedy)
    outputPath = f"{dirPath}_txt/{n:05d}.txt"
    with open(outputPath, "w") as file:
        for x in range(dArrayx.shape[0]):
            for y in range(dArrayx.shape[1]):
                file.write(f"{x} {y} {dArrayx[x, y]} {dArrayy[x, y]}\n")
    file.close()
    n += 1