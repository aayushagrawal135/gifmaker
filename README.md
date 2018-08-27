Navcam images from MCS are added with contrast, bundled together as frames, and put as a gif. 


# Processing Navcam images to output a GIF

[dataset used](https://mars.nasa.gov/msl/multimedia/raw/?s=1138&camera=NAV%5FRIGHT%5F) Aim is to enhanced or lighten the image. 

## Getting Started


### Prerequisites

Libraries to be installed: numpy, imageio, cv2
[install for cv2](https://medium.com/@debugvn/installing-opencv-3-3-0-on-ubuntu-16-04-lts-7db376f93961)


## Running the files

From the progs folder, keep all the files in it in the same directory. A sample image as <van.png> is present in the <progs> folder. Set the input file name according in eishi.py. Run it as <python3 eishi.py> . Python 3.5.2 used. 



### support.py file

Processing functions namely <findMed, vSize, zStretch, toGray, meanSubtracion> are implemented.

findMed: It makes buckets recursively find finding median. And now left and and right bound for the bucket are found by taking median in the left half and the right half. 

vSize: Reshapes a (m, n) to (1, m*n).

zStretch: Main handler for findMed.

toGray: Reduce 3 channel image to 1 chaannel image.

meanSubtraction: Normalizes the pixel values. 
