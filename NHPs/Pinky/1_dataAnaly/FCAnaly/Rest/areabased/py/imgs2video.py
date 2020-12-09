import os
import cv2
from tkinter import Tk
from tkinter import filedialog


def generate_video(genvideofile = None, images = None): 
    """

        args:
            genvideofile: the generated video file full path

            images: the required images

    
    """

    if images is None:
        Tk().withdraw()
        images = filedialog.askopenfilenames(title='choose images')


    if genvideofile is None:
        imagesfolder = os.path.split(images[0])[0]
        genvideofile = os.path.join(imagesfolder, 'combined.avi')

    frame = cv2.imread( images[0])
  
    # setting the frame width, height width 
    # the width, height of first image 
    height, width, layers = frame.shape   
  
    video = cv2.VideoWriter(genvideofile, 0, 1, (width, height))  
  
    # Appending the images to the video one by one 
    for image in images:  
        video.write(cv2.imread(image))  
      
    # Deallocating memories taken for window creation 
    cv2.destroyAllWindows()  
    video.release()  # releasing the video generated 
    print("Generated video")
            

if __name__ == "__main__":
    generate_video(genvideofile = None, images = None)