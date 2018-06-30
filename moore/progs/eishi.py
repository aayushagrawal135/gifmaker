import support as sup
import cv2
import imageio

def makegif(filenames):
	gifname = "/home/aayush/Downloads/Curiosity/moore/set1/movie.gif"
	images = []
	for file in filenames:
		images.append(file)
	imageio.mimsave(gifname, images)

filename = "/home/aayush/Downloads/Curiosity/moore/set1/"
images = []
for x in range(8):
	name = filename + str(x) + ".jpg"
	img = cv2.imread(name)
	img = sup.toGray(img)
	img = sup.zStretch(img, 3)
	images.append(img)
makegif(images)