import sys
import struct
import numpy as np
import scipy.misc
from tifffile import imsave
import numpy as np

inFile = sys.argv[1]

f=open(inFile,"rb")
(num,) = struct.unpack('i', f.read(4))
dimX = abs(num)+1
(num,) = struct.unpack('i', f.read(4))
dimY = abs(num)+1
(num,) = struct.unpack('i', f.read(4))
dimZ = abs(num)+1

(num,) = struct.unpack('d', f.read(8))
(num,) = struct.unpack('d', f.read(8))
(num,) = struct.unpack('d', f.read(8))
(num,) = struct.unpack('d', f.read(8))
(num,) = struct.unpack('d', f.read(8))
(num,) = struct.unpack('d', f.read(8))

vol = np.zeros((dimX, dimY, dimZ))
for k in range(0, dimZ):
	for j in range(0, dimY):
		for i in range(0, dimX):
			(num,) = struct.unpack('f', f.read(4))
			vol[i,j,k] = -num;

f.close()

imsave(sys.argv[2], vol.astype(np.float32))
