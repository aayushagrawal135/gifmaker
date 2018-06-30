import support as sup
import cv2
import imageio

def makegif(filenames):
	gifname = "/home/aayush/Downloads/Curiosity/moore/set1/movie.gif"
	images = []
	for file in filenames:
		images.append(file)
	imageio.mimsave(gifname, images, fps = 3)

#---------------------------------------------------------------

filename = "/home/aayush/Downloads/Curiosity/moore/set1/"
images = []

for i in range(8):
	name = filename + str(i) + ".jpg"
	img = cv2.imread(name)
	img = sup.toGray(img)
	images.append(img)

images = sup.meanSubtraction(images)
res = []
for i in range(8):
	img = sup.zStretch(images[i], 3)
	res.append(img)

makegif(res)
