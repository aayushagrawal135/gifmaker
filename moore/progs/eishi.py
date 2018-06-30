import support as sup
import cv2
import imageio

def makegif(filenames, gifname):
	images = []
	for file in filenames:
		images.append(file)
	imageio.mimsave(gifname, images, fps = 3)

#---------------------------------------------------------------

def wrapper(foldername, gifname, iter = 3, fps = 3):
	images = []

	for i in range(8):
		name = foldername + str(i) + ".jpg"
		img = cv2.imread(name)
		img = sup.toGray(img)
		images.append(img)

	images = sup.meanSubtraction(images)
	res = []
	for i in range(8):
		img = sup.zStretch(images[i], 3)
		res.append(img)

	makegif(res, gifname)
#-------------------------------------------------

foldername = "/home/aayush/Downloads/Curiosity/moore/set1/"
gifname = "/home/aayush/Downloads/Curiosity/moore/set1/movie.gif"

#Assumed: There are atleast 8 images names as x.jpg where x lies from 0 to number of images. 
wrapper(foldername, gifname)