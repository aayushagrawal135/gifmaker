import cv2
import numpy as np
import os

filename = "NRB_498502911EDR_M0500676NCAM00536M_-thm.jpg"
img = cv2.imread(filename)

print(os.stat(filename).st_size)
