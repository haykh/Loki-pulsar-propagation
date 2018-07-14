import struct
import numpy as np

def loki_read(fname, ncol = 4, d_type = 8):
    """
    Function to read binary data by `Loki` pulsar code, analize and return the results.
    Typically reads 8 byte (`d_type`) double.

    Parameters
    ----------
    fname : string
        full path to the binary file
    ncol : integer (4 by default)
        number of columns to read
    d_type : integer (8 by default)
        size in # of bytes for the data in binary
        number of bytes: double = 8, int = 4, float = 4

    Returns
    -------
    array_like
        array of floats containing `ncol` columns

    Raises
    ------
    ArithmeticError
        wrong number of columns, `ncol` for a given data in `fname`
    """
    with open(fname, mode='rb') as file: # b is important -> binary
        fileContent = file.read()
        data = struct.unpack("d" * (len(fileContent) // 8), fileContent) # 8 -> sizeof(double)
    if (len(data) % ncol != 0):
        raise ArithmeticError("ERROR: Enter the correct ncol")
    data = np.reshape(np.array(data), (int(len(data) / ncol), ncol)) # ncol -> number of columns
    return data[data[:,0].argsort()]

def fix_PA(PA_data, threshold = 300):
    """
    Shifts P.A. data [in degrees] by 360 if necessary

    Parameters
    ----------
    PA_data : array_like
    threshold : float
        threshold in degrees defining how much of a difference is worth shifting down

    Returns
    -------
    array_like
        shifted P.A. data array
    """
    new_PA_data = np.copy(PA_data)
    for i in range(1, len(new_PA_data)):
        if np.abs(new_PA_data[i - 1] - new_PA_data[i]) > threshold:
            new_PA_data[i] -= 360
    return new_PA_data
