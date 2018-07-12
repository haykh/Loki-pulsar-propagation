import struct
import numpy as np

def loki_read(fname, ncol = 4):
    with open(fname, mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        data = struct.unpack("d" * (len(fileContent) // 8), fileContent) # 8 -> sizeof(double)
    data = np.reshape(data, (len(data) / ncol, ncol)) # ncol -> number of columns
    return data[data[:,0].argsort()]
