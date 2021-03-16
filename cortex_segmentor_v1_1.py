#!/opt/local/bin/python

__author__ = "Kai Dierkes"
__date__ = "2017"
__maintainer__ = "Andrew G. Clark"
__email__ = "andrew.clark@curie.fr"

"""

This GUI-based program is designed to read images, perform segmentations
of an enriched region in the cell periphery and generate linescans normal
to the segmentation contour.

Images and segmentation files can be loaded, and segmentation parameters
can be saved. Linescans can also be generated and are automatically saved
in the same directory as the image.

Classes:
    Experiment
    App

"""

import os
import math

import numpy as np
from scipy.optimize import leastsq
from tkinter import *
from tkinter.filedialog import *
master = Tk()
from PIL import Image, ImageTk
import matplotlib.pyplot as plt

def fitfunc(p,phi):
    """Fourier series to fit the segmentation contour to the selected fitpoints."""

    value = p[0]

    for i in range(1,len(p),2):

        k = (i+1)/2
        value += p[i]*np.sin(k*phi)+p[i+1]*np.cos(k*phi)

    return value

def errfunc(p,phi,xdat):
    """Error function used for fitting"""

    return fitfunc(p,phi)-xdat

class Experiment:
    """A class to contain the image and segmentation data
    for the current image.

    """

    def __init__(self, file):
        """Initializes the experiment class

        Args:
            file (str): file path
        """

        self.filename = file
        self.segmentation = ""
        self.directory = os.path.dirname(file)
        self.cellimage  = Image.open(self.filename)
        self.cellimagematrix = self.cellimage.load()
        self.size_x = self.cellimage.size[0]
        self.size_y = self.cellimage.size[1]
        self.no_pixels = self.size_x*self.size_y
        self.no_stacks = self.countstacks(self.cellimage)

        self.cellimage.seek(0)
        self.current_layer = 0

        self.fit_x = []
        self.fit_y = []
        self.fit_points = []
        self.outline_pixel_list = []
        self.linescan_pars = []

        for i in range(0,self.no_stacks,1):

            self.fit_x.append([False,[]])
            self.fit_y.append([False,[]])
            self.fit_points.append([])
            self.outline_pixel_list.append([])
            self.linescan_pars.append([0,360,50,50])

    def countstacks(self,image):
        """Counts the number of slices in the image.

        Args:
            image (PIL Image): image to be counted

        Returns:
            stack_counter (int): number of slices in the image stack

        """
         
        stack_counter = 0
        eof_indicator = 0
         
        while eof_indicator!=1:
            try:
                image.seek(stack_counter)
                stack_counter += 1
            except EOFError:
                  eof_indicator = 1
         
        return stack_counter

    def change_fit_points(self,xindex,yindex,fit_toggle,radius=2):
        """Adds or removes fit points from the image

        Args:
            xindex (int): x-coordinate of location of fit point addition/removal
            yindex (int): y-coordinate of location of fit point addition/removal
            fit_toggle (bool): toggle to add (==1) or remove (==0) fit point

        """
         
        for i in range(-radius,radius+1,1):
            for j in range(-radius,radius+1,1):
                if  0<=xindex+i<self.size_x and 0<=yindex+j<self.size_y:

                    if fit_toggle==1:
                        if [xindex+i,yindex+j] not in self.fit_points[self.current_layer]:
                            self.fit_points[self.current_layer].append([xindex+i,yindex+j])

                    else:
                        while [xindex+i,yindex+j] in self.fit_points[self.current_layer]:
                            self.fit_points[self.current_layer].remove([xindex+i,yindex+j])

    def seek(self,layer):
        """Changes current layer (slice) of the image

        Args:
            layer (int): desired layer to make current

        """
         
        if layer>=0 and layer<self.no_stacks:
            self.cellimage.seek(layer)
            self.cellimagematrix = self.cellimage.load()
            self.current_layer = layer
            
        else:
            return 0

class App:
    """A class for the GUI"""

    def __init__(self,root):
        """Initializes the GUI App Class with windows, buttons, etc.

        Args:
            root (Tkinter.Tk instance): Tk instance for drawing GUI
        """

        self.frame = Frame(root,width=60, height=512)
        self.frame.grid(row=0,column=0,padx=10)                  

        self.openbutton = Button(self.frame, text="Open file",command = self.openfile,width=15)
        self.openbutton.grid(row=0,column=0)
        
        self.loadbutton = Button(self.frame, text="Load", command=self.load,width=15)
        self.loadbutton.grid(row=0,column=1)
        
        self.seekplus = Button(self.frame, text="Next layer", command = lambda : self.seek(self.cell.current_layer+1),width=15)
        self.seekplus.grid(row=10,column=1)
        
        self.seekminus = Button(self.frame, text="Previous layer",command = lambda : self.seek(self.cell.current_layer-1),width=15)
        self.seekminus.grid(row=10,column=0)    
            
        self.image_toggle = IntVar()
        self.fit_points_toggle = IntVar()
        self.segmentation_toggle = IntVar()
        self.linescan_toggle = IntVar()
        self.linescanimage_toggle = IntVar()
        
        self.toggles = [[self.image_toggle,"Image"],
                        [self.fit_points_toggle,"Fit points"],
                        [self.segmentation_toggle,"Segmentation"],
                        [self.linescan_toggle,"Linescan - Points"],
                        [self.linescanimage_toggle,"Linescan - Image"]]

        for i in range(0,len(self.toggles),1):
                   
            self.toggles[i].append(Checkbutton(self.frame, variable = self.toggles[i][0],command = self.draw))

            if i!=4:
                self.toggles[i][2].select()
            else:
                self.toggles[i][2].deselect()

            self.toggles[i][2].grid(row=i+20,column=1)
            self.toggles[i].append(Label(self.frame,text=self.toggles[i][1]))
            self.toggles[i][3].grid(row=i+20,column=0)
    
        self.phistart = Entry(self.frame,bg="white",width=10)
        self.phistart.insert(0,"0")
        self.phistart.grid(row=50,column=1)
    
        self.phiend = Entry(self.frame,bg="white",width=10)
        self.phiend.insert(0,"360")
        self.phiend.grid(row=51,column=1)
    
        self.outerlength = Entry(self.frame,bg="white",width=10)
        self.outerlength.insert(0,"50")
        self.outerlength.grid(row=52,column=1)    
    
        self.innerlength = Entry(self.frame,bg="white",width=10)
        self.innerlength.insert(0,"50")
        self.innerlength.grid(row=53,column=1)
        
        self.resetradius = Entry(self.frame,bg="white",width=10)
        self.resetradius .insert(0,"2")
        self.resetradius.grid(row=65,column=0)
        
        self.resetlabel = Label(self.frame,text="Fitpoint radius")
        self.resetlabel.grid(row=65,column=1)
        
        self.line_scan_labels = [["Startangle"],["Endangle"],["Outer length"],["Inner length"]]
        
        for i in range(0,len(self.line_scan_labels),1):
           
           self.line_scan_labels[i].append(Label(self.frame,text=self.line_scan_labels[i][0]))
           self.line_scan_labels[i][1].grid(row=50+i,column=0)
        
        self.copylinescanbutton = Button(self.frame, text="Copy to all",command = self.copy_linescan_parameters,width=6)
        self.copylinescanbutton.grid(row = 50,column=2)
        
        self.linescanbutton = Button(self.frame, text="Linescan",command = self.linescan,width=15)
        self.linescanbutton.grid(row = 54)
        
        self.linescanshowbutton = Button(self.frame, text="Show endpoints",command = self.linescan_show_endpoints,width=15)
        self.linescanshowbutton.grid(row = 54,column=1)
      
        self.scale_fit_points = Scale(self.frame, from_=0, to=0.1, orient=HORIZONTAL,resolution=0.0005,troughcolor="white",showvalue=1,borderwidth=2)
        self.scale_fit_points.set(0.075)
        self.scale_fit_points.grid(row=30,column=1)
        
        self.fixpointsbutton = Button(self.frame, text="Choose fit points",command =self.choose_fit_points,width=15)
        self.fixpointsbutton.grid(row=30,column=0)
        
        self.segmentbutton = Button(self.frame, text="Segment",command = lambda: self.segment(-1),width=15)
        self.segmentbutton.grid(row=40)
        
        #self.segmentmodeslabel = Label(self.frame,text = "No. of modes")
        #self.segmentmodeslabel.grid(row=40,column = 0)
        
        self.segmentmodes = Entry(self.frame, bg="white",width=10)
        self.segmentmodes.insert(0,"10")
        self.segmentmodes.grid(row=40,column=1)
        
        self.hr1 = Canvas(self.frame,height=3,width=200,bd=-2)
        self.hr1.grid(row=43,columnspan=2,pady =2)
        
        self.copysegmenationentry = Entry(self.frame,bg="white",width=10)
        self.copysegmenationentry.insert(0,"10")
        self.copysegmenationentry.grid(row=45,column=1)
        
        self.segmentcopybutton = Button(self.frame, text="Copy segmentation",command = self.copysegmentation,width=15)
        self.segmentcopybutton.grid(row=45,column=0)    
        
        self.segmentcopychoicewhat = StringVar(root)
        self.segmentcopychoicewhat.set("Entry")                   
        
        self.segmentcopychoice = OptionMenu(self.frame,self.segmentcopychoicewhat,"Entry", "Odd+1->Even","Even-1->Odd", "Ch1->Ch2/3", "Ch2->Ch1/3", "Ch3->Ch1/2")
        self.segmentcopychoice.grid(row=45,column=2)
        
        self.analyzewhatentry =  Entry(self.frame,bg="white",width=10)
        #self.analyzewhatentry.insert(0,"")
        self.analyzewhatentry.grid(row=26,column=1)
        
        self.analyzewhatlabel = Label(self.frame, text="Mode:",padx=5,pady = 5)
        self.analyzewhatlabel.grid(row=26,column=0)
        
        self.analyzewhat = StringVar(root)
        self.analyzewhat.set("Current") 

        self.analyzewhatchoice = OptionMenu(self.frame, self.analyzewhat, "Current","Entry","All","Even","Odd",
                                            "Every third (from 0)", "Every third (from 1)", "Every third (from 2)",
                                            "All from entry","Every second from entry")
        self.analyzewhatchoice.grid(row=26,column=2)

        self.hr1 = Canvas(self.frame,height=3,width=200,bd=-2)
        self.hr1.grid(row=25,columnspan=2,pady =2)
        
        self.resetfitpointsbutton = Button(self.frame, text="Reset Fitpoints", relief="groove",command=self.resetfitpoints,width=15)
        self.resetfitpointsbutton.grid(row=65,column=2)
        
        self.button = Button(self.frame, text="Quit", relief="groove",command=self.frame.quit,width=15)
        self.button.grid(row=70,column=0)
        
        self.savebutton = Button(self.frame, text="Save", command=self.save,width=15)
        self.savebutton.grid(row=70,column=1)
            
        
        for i in range(10,80,10):
        
            self.hr1 = Canvas(self.frame,height=3,width=200,bd=-2)
            self.hr1.grid(row=i-1,columnspan=2,pady =2)
                  
                
        self.display = Canvas(root, width=512, height=512, background="black",bd=-1)       
        self.display.bind("<Button-1>", self.add_fit_points)
        self.display.bind("<Button-2>", self.remove_fit_points)
        self.display.grid(row=0,column=2)
        
        self.cell = Experiment("./startupscreen.tif")

        self.image_canvas = ImageTk.PhotoImage(self.cell.cellimage)
        self.display.create_image((self.cell.size_x/2,self.cell.size_y/2),image=self.image_canvas) 
        
        self.frame.bind_all("<Down>", self.down)
        self.frame.bind_all("<Left>", self.left)
        self.frame.bind_all("<Right>",self.right)

    def resetfitpoints(self):
        """Resets (deletes) all current fit points"""
        
        self.cell.fit_points[self.cell.current_layer] = []
        self.draw()  
        self.display.update_idletasks()
        
    def down(self,event):
        """Segments the current cell layer.

        Args:
            event (Tk keystroke event)

        """
        self.segment(self.cell.current_layer)
    
    def left(self,event):
        """Seeks one slice backward

        Args:
            event (Tk keystroke event)

        """
        self.seek(self.cell.current_layer-1)
    
    def right(self,event):
        """Seeks one slice forward

        Args:
            event (Tk keystroke event)

        """
        self.seek(self.cell.current_layer+1)
           
    def copy_linescan_parameters(self):
        """Copies linescan parameters from one frame to another"""
       
        outer = float(self.outerlength.get())
        inner = float(self.innerlength.get())
        phistart = float(self.phistart.get())
        phiend = float(self.phiend.get())
         
        for layer in range(0,self.cell.no_stacks):
            self.cell.linescan_pars[layer] = [phistart, phiend, outer, inner]
    
    def seek(self,layer):
        """Seeks image and segmentation/linescan data to a specific layer (slice)

        Args:
            layer (int): desired slice

        """
    
        self.store_linescan_entries()
        self.cell.seek(layer)
        self.load_linescan_entries()
        self.draw()
      
    def generate_framelist(self):
        """Makes a list of frames to process based on users' selection in the GUI

        Returns:
            layers (list): list of layers to be processed
        """
           
        analyze_choice = self.analyzewhat.get()
        analyze_list   = self.analyzewhatentry.get()

        if "'" in analyze_list:
            liste = analyze_list.split(",")

        else:
            liste = analyze_list.split()
                  
        if len(liste)==0:
            liste = [0]
                     
        if analyze_choice=="Current":
            layers = [self.cell.current_layer]

        elif analyze_choice=="Entry":
            layers = []
            for entry in liste:
                layers.append(int(entry))

        elif analyze_choice=="All":
            layers = range(0,self.cell.no_stacks,1)

        elif analyze_choice=="Even":
            layers = range(0,self.cell.no_stacks,2)

        elif analyze_choice=="Odd":
            layers = range(1,self.cell.no_stacks,2)

        elif analyze_choice=="Every third (from 0)":
            layers = range(0,self.cell.no_stacks,3)

        elif analyze_choice=="Every third (from 1)":
            layers = range(1,self.cell.no_stacks,3)

        elif analyze_choice=="Every third (from 2)":
            layers = range(2,self.cell.no_stacks,3)

        elif analyze_choice=="All from entry":
            if 0<=int(liste[0]) and int(liste[0])<self.cell.no_stacks:
                layers = range(int(liste[0]),self.cell.no_stacks,1)
            else:
                layers = [0]

        elif analyze_choice=="Every second from entry":
            if 0<=int(liste[0]) and int(liste[0])<self.cell.no_stacks:
                layers = range(int(liste[0]),self.cell.no_stacks,2)
            else:
                layers = [0]
                      
        print("Layers to work on:",layers)
        return layers
            
    def copysegmentation(self):
        """Copies segmentation and linescan parameters from one frame to another"""

        #prepares lists for odds/evens and ch1-3
        odds  = range(1,self.cell.no_stacks,2)
        evens = range(0,self.cell.no_stacks,2)

        ch1 = range(0,self.cell.no_stacks,3)
        ch2 = range(1,self.cell.no_stacks,3)
        ch3 = range(2,self.cell.no_stacks,3)

        #gets copy choice mode and determines list of copies
        mode = self.segmentcopychoicewhat.get()
        listofcopies = []

        if mode=="Entry":
            noframe = int(self.copysegmenationentry.get())
            listofcopies = [[noframe,self.cell.current_layer]]

        elif mode=="Odd+1->Even":
            for i in range(0,len(odds),1):
                listofcopies.append([odds[i],evens[i]])

        elif mode=="Even-1->Odd":
            for i in range(0,len(odds),1):
                listofcopies.append([evens[i],odds[i]])

        elif mode=="Ch1->Ch2/3":
            for i in range(0,len(ch1),1):
                listofcopies.append([ch1[i],ch2[i]])
                listofcopies.append([ch1[i],ch3[i]])

        elif mode=="Ch2->Ch1/3":
            for i in range(0,len(ch1),1):
                listofcopies.append([ch2[i],ch1[i]])
                listofcopies.append([ch2[i],ch3[i]])

        elif mode=="Ch3->Ch1/2":
            for i in range(0,len(ch1),1):
                listofcopies.append([ch3[i],ch1[i]])
                listofcopies.append([ch3[i],ch2[i]])

        else:
            raise ValueError("Not a valid copy mode chioce!")
            
        self.store_linescan_entries()

        #copies fit points, segmentation and linescan parameters
        for movement in listofcopies:
            
            toframe = movement[1]
            noframe = movement[0]
            
            self.cell.fit_x[toframe][0]=self.cell.fit_x[noframe][0]
            self.cell.fit_x[toframe][1] = []
               
            for entry in  self.cell.fit_x[noframe][1]:
                  self.cell.fit_x[toframe][1].append(entry)
               
            self.cell.fit_y[toframe][0] = self.cell.fit_y[noframe][0]
            self.cell.fit_y[toframe][1] = []
               
            for entry in  self.cell.fit_y[noframe][1]:
                  self.cell.fit_y[toframe][1].append(entry)
               
            self.cell.fit_points[toframe] = []
               
            for entry in self.cell.fit_points[noframe]:
               self.cell.fit_points[toframe].append(entry)
               
            self.cell.linescan_pars[toframe] = []
               
            for entry in self.cell.linescan_pars[noframe]:
                self.cell.linescan_pars[toframe].append(entry)
         
        self.load_linescan_entries()
        self.draw()
                    
    def save(self):
        """Saves current fit points, segmentation and linescan parameters for the whole image stack"""
        
        self.store_linescan_entries()
        current = os.getcwd()
        os.chdir(self.cell.directory)
                 
        output = asksaveasfile(mode='w',filetypes=[("csd", "*.csd")],initialfile="%s"%self.cell.segmentation)

        #goes through slices, collects and writes the data
        for i in range(0,self.cell.no_stacks,1):

            if self.cell.fit_x[i][0]:
                output.write("1\n")
            else:
                output.write("0\n")
                  
            for k in range(0,len(self.cell.fit_x[i][1]),1):
                output.write("%e "%self.cell.fit_x[i][1][k])
                   
            output.write("\n")
             
            for k in range(0,len(self.cell.fit_x[i][1]),1):
                output.write("%e "%self.cell.fit_y[i][1][k])
             
            output.write("\n")
                
            for k in range(0,4,1):
                output.write("%i "%self.cell.linescan_pars[i][k])
             
            output.write("\n")
             
            for point in self.cell.fit_points[i]:
                output.write("%i %i "%(point[0],point[1]))
                  
            output.write("\n")
             
        output.close()
        os.chdir(current)
                  
    def load(self):
        """Loads a new segmentation file"""
    
        data_file =  askopenfilename(filetypes=[("csd", "*.csd")],initialfile="%s"%self.cell.segmentation)
         
        try:
            ifile = open(data_file,"r")
            self.cell.segmentation = data_file.split("/")[-1]
        except:
            return
         
        dat_temp = ifile.readlines()
        ifile.close()

        #adds the segmentation data to the currently open image
        for i in range(0,self.cell.no_stacks,1):
                        
            if int(dat_temp[i*5+0].split()[0])==1:
               self.cell.fit_x[i][0] = True
               self.cell.fit_y[i][0] = True
               
            else:
               self.cell.fit_x[i][0] = False
               self.cell.fit_y[i][0] = False
               
            self.cell.fit_x[i][1] = []
               
            for fourier in dat_temp[i*5+1].split():
                self.cell.fit_x[i][1].append(float(fourier))
               
            self.cell.fit_y[i][1] = []
               
            for fourier in dat_temp[i*5+2].split():
                self.cell.fit_y[i][1].append(float(fourier))
            
            linescan_pars = dat_temp[i*5+3].split()
                        
            for j in range(0,4,1):
                self.cell.linescan_pars[i][j]=int(linescan_pars[j])
            
            self.cell.fit_points[i] = []
            
            fit_points_list = dat_temp[i*5+4].split()
            
            for j in range(0,len(fit_points_list),2):
                self.cell.fit_points[i].append([int(fit_points_list[j]),int(fit_points_list[j+1])])
               
            self.cell.outline_pixel_list[i] = []   
         
        self.load_linescan_entries()
        self.draw()
       
    def load_linescan_entries(self):
        """Loads the linescan parameters to the GUI"""
          
        outer = self.cell.linescan_pars[self.cell.current_layer][2]
        inner = self.cell.linescan_pars[self.cell.current_layer][3]
        phistart = self.cell.linescan_pars[self.cell.current_layer][0]
        phiend = self.cell.linescan_pars[self.cell.current_layer][1]
         
        self.outerlength.delete(0,END)
        self.innerlength.delete(0,END)
        self.phistart.delete(0,END)
        self.phiend.delete(0,END)

        self.outerlength.insert(0,"%i"%outer)
        self.innerlength.insert(0,"%i"%inner)
        self.phistart.insert(0,"%i"%phistart)
        self.phiend.insert(0,"%i"%phiend)

        self.draw()
    
    def store_linescan_entries(self):
        """Gets the linescan parameters from the GUI"""
  
        outer = float(self.outerlength.get())
        inner = float(self.innerlength.get())
        phistart = float(self.phistart.get())
        phiend = float(self.phiend.get())

        self.cell.linescan_pars[self.cell.current_layer] = [phistart,phiend,outer,inner]
    
    def choose_fit_points(self):
        """Automatically selects fit points for fitting the segmentation contour
        based on high thresholding and sequential deletion/refinement of the fit points.

        The quality of the automatic segmentation can be adjusted to fit needs for specific
        images here by changing the 'algorithm' list. This list determines what points should
        be removed at each iteration. Each iteration is a sublist of three values indicating:
        [direction (-1 = outside the current fit, 1 = inside the current fit),
         number of modes used for the fit after this iteration,
         the minimum distance away from the segmentation contour for fit point deletion]
        """
       
        layers = self.generate_framelist()

        #loops through the layers to be processed and performs the automated fit point selection
        for layer in layers:
             
            self.cell.seek(layer)
            
            current_modes  = self.segmentmodes.get()
            self.segmentmodes.delete(0,END)
            self.segmentmodes.insert(0,"4")
            
            self.cell.fit_points[self.cell.current_layer] = []
            
            #get a sorted list of pixel intensities
            intensities  = []
            for i in range(0,self.cell.size_x,1):
                for k in range(0,self.cell.size_y,1):
                    intensities.append([np.sum(self.cell.cellimagematrix[i,k]),i,k])
            intensities.sort()

            #initially chooses fit points based on simple threshold (factor selected in GUI)
            for i in range(1,int(self.scale_fit_points.get()*self.cell.no_pixels),1):
                self.cell.fit_points[self.cell.current_layer].append([intensities[-i][1],intensities[-i][2]])
            
            self.draw()  
            self.display.update_idletasks()

            #sequentially and progressively deletes fit points that are far away from the fit
            algorithm = [[-1,4,25],
                         [1,4,25],
                         [-1,4,10],
                         [1,4,15],
                         [-1,8,10],
                         [1,8,15],
                         [-1,8,3],
                         [-1,10,2]]

            #other examples of algorithms to use for automated fit point selection
            # algorithm = [[-1,4,5]
            #              [1,4,25],
            #              [-1,4,25],
            #              [-1,4,10],
            #              [1,4,10],
            #              [1,8,10],
            #              [-1,8,10],
            #              [1,8,4],
            #              [1,8,2],
            #              [-1,8,4],
            #              [1,10,4]]
            # algorithm = [[ 1,4,25],
            #              [1,4,25],
            #              [-1,4,45],
            #              [1,4,15],
            #              [1,4,10],
            #              [-1,6,15],
            #              [1,6,10],
            #              [1,6,8],
            #              [-1,8,15],
            #              [-1,8,15]]
            # algorithm = [[1,8,10],
            #              [1,8,7],
            #              [-1,8,7],
            #              [-1,10,12],
            #              [1,10,8],
            #              [-1,10,10],
            #              [1,10,10],
            #              [-1,11,9]]

            #iterates through the algorithm list and removes fit points
            for iteration in algorithm:
               
                side = iteration[0]
                no_modes = iteration[1]
                distance = iteration[2]

                self.segmentmodes.delete(0,END)
                self.segmentmodes.insert(0,"%s"%no_modes)

                self.segment(self.cell.current_layer)

                delete_list = []

                for entry in self.cell.fit_points[self.cell.current_layer]:
                  
                    PHI = math.atan2(float(256-entry[1]),float(256-entry[0]))
                    i = int(round(fitfunc(self.cell.fit_x[self.cell.current_layer][1],PHI)))
                    j = int(round(fitfunc(self.cell.fit_y[self.cell.current_layer][1],PHI)))
                  
                    if math.sqrt((i-entry[0])**2+(j-entry[1])**2)>distance and 0<side*(math.sqrt((256-i)**2+(256-j)**2)-math.sqrt((256-entry[0])**2+(256-entry[1])**2)):
                     
                        delete_list.append(entry)
               
                for entry in delete_list:
                  
                    self.cell.fit_points[self.cell.current_layer].remove(entry)

            #re-segments and puts modes back to what it was before
            self.segment(self.cell.current_layer)
            self.segmentmodes.delete(0,END)
            self.segmentmodes.insert(0,"%s"%current_modes)
            self.segment(self.cell.current_layer)

            self.draw()
            self.display.update_idletasks()

    def draw(self): 
        """Draws the cell image, fit points and segmentation (if toggled in GUI)"""

        allitems = self.display.find_all()
          
        for item in allitems:
            self.display.delete(item)

        #displays image
        if self.image_toggle.get() == 1:
            self.image_canvas = ImageTk.PhotoImage(self.cell.cellimage)
            self.display.create_image((self.cell.size_x/2,self.cell.size_y/2),image=self.image_canvas)

        #displays fit points
        if self.fit_points_toggle.get() == 1:
            for point in self.cell.fit_points[self.cell.current_layer]:
                self.display.create_rectangle((point[0],point[1],point[0],point[1]),width=0,fill="green")

        #displays segmentation
        if self.segmentation_toggle.get() == 1:
            if self.cell.fit_x[self.cell.current_layer][0]:
                                         
                cell_outline = []

                for PHI in np.arange(0,2*np.pi,0.01):
                    i = int(round(fitfunc(self.cell.fit_x[self.cell.current_layer][1],PHI)))
                    j = int(round(fitfunc(self.cell.fit_y[self.cell.current_layer][1],PHI)))
                    cell_outline.append([i,j])
      
                self.display.create_line(cell_outline,fill="red",width="2.0")
                                          
        if self.linescan_toggle.get() == 1:
                                
            draw_list = []

            for entry in self.cell.outline_pixel_list[self.cell.current_layer]:
                draw_list.append([entry[0],entry[1]])
    
            try:
                self.display.create_line(draw_list,fill="white",width="4.0")
            except:
                pass
         
        self.display.create_text((20,20),text = "Layer %i/%i"%(self.cell.current_layer,self.cell.no_stacks),fill="white",anchor="w")
        self.display.create_text((256,500),text = "%s"%self.cell.filename.split("/")[-1],fill="white")
        self.display.create_text((256,480),text = "%s"%self.cell.segmentation,fill="white")
        #self.display()
        self.display.update_idletasks()

    def openfile(self):
        """Opens image file"""

        image_file = askopenfilename(filetypes=[("tif", "*.tif")], initialdir=self.cell.directory)
        self.cell = Experiment(image_file)
      
        outer = self.cell.linescan_pars[self.cell.current_layer][2]
        inner = self.cell.linescan_pars[self.cell.current_layer][3]
        phistart = self.cell.linescan_pars[self.cell.current_layer][0]
        phiend = self.cell.linescan_pars[self.cell.current_layer][1]

        self.outerlength.delete(0,END)
        self.innerlength.delete(0,END)
        self.phistart.delete(0,END)
        self.phiend.delete(0,END)

        self.outerlength.insert(0,"%i"%outer)
        self.innerlength.insert(0,"%i"%inner)
        self.phistart.insert(0,"%i"%phistart)
        self.phiend.insert(0,"%i"%phiend)

        self.draw()
                
    def add_fit_points(self,event):
       """Adds fit point to current mouse location

       Args:
           event (Tk mouse click event)

       """
       radius = int(self.resetradius.get())
       self.cell.change_fit_points(event.x,event.y,1,radius)
       self.draw()
       #print('clicked to add at', event.x, event.y)

    def remove_fit_points(self,event):
       """Adds fit point to current mouse location

       Args:
           event (Tk mouse click event)

       """
       radius = int(self.resetradius.get())
       self.cell.change_fit_points(event.x,event.y,0,radius)
       self.draw()
       #print('clicked to remove at', event.x, event.y)

    def segment(self,mode):
        """Fits and draws a segmentation contour based on current fit points

        Args:
            mode (int): instruction to generate frame list for processing multiple frames
            if mode==-1, otherwise, mode corresponds to the single frame to be processed

        """
        if mode==-1:
            scanlayers = self.generate_framelist()
        else:
            scanlayers = [mode]
     
        for layer in scanlayers:
                  
            self.cell.seek(layer)
            self.draw()
            self.display.update_idletasks()
            
            no_fit_points = len(self.cell.fit_points[self.cell.current_layer])
            
            #gets angles of fit points
            phi = [] #zeros(no_fit_points)
            x   = [] #zeros(no_fit_points)
            y   = [] #zeros(no_fit_points)
                     
            for i in range(0,no_fit_points,1):
                           
                xi = self.cell.fit_points[self.cell.current_layer][i][0]
                yi = self.cell.fit_points[self.cell.current_layer][i][1]

                i0 = self.cell.fit_points[self.cell.current_layer][i][0]
                j0 = self.cell.fit_points[self.cell.current_layer][i][1]

                #for j in range(0,self.cell.cellimage.getpixel((i0,j0)),1):

                phi.append(math.atan2(float(256-yi),float(256-xi)))
                x.append(xi)
                y.append(yi)

            phi = np.reshape(phi,len(phi))
            x = np.reshape(x,len(x))
            y = np.reshape(y,len(y))

            #performs fit
            no_modes = 2*(int(self.segmentmodes.get())-1)+1
            p0 = np.ones(no_modes)
                  
            self.cell.fit_x[self.cell.current_layer][1],success  = leastsq(errfunc, p0, args=(phi,x), maxfev = 500)
            self.cell.fit_x[self.cell.current_layer][0] = True
            self.cell.fit_y[self.cell.current_layer][1],success  = leastsq(errfunc, p0, args=(phi,y), maxfev = 500)
            self.cell.fit_y[self.cell.current_layer][0] = True                                                                                                   

            self.draw()
            self.display.update_idletasks()
        
    def linescan_show_endpoints(self):
        """Displays the endpoints of the linescan as determined by start and end
        angles in the GUI

        """
        phistart = float(self.phistart.get())
        phiend = float(self.phiend.get())

        istart = int(np.floor(fitfunc(self.cell.fit_x[self.cell.current_layer][1],phistart*np.pi/180.0)))
        jstart = int(np.floor(fitfunc(self.cell.fit_y[self.cell.current_layer][1],phistart*np.pi/180.0)))

        normalistart, normaljstart = self.normal_vector(phistart*np.pi/180.0)

        iend = int(np.floor(fitfunc(self.cell.fit_x[self.cell.current_layer][1],phiend*np.pi/180.0)))
        jend = int(np.floor(fitfunc(self.cell.fit_y[self.cell.current_layer][1],phiend*np.pi/180.0)))

        normaliend, normaljend = self.normal_vector(phiend*np.pi/180.0)

        self.draw()

        #length = 30

        il = float(self.innerlength.get())
        ol = float(self.outerlength.get())

        #draws line
        # length = 30
        # self.display.create_line([istart+length*normalistart,jstart+length*normaljstart,istart-length*normalistart,jstart-length*normaljstart],width=2,fill="grey")
        # self.display.create_line([iend+length*normaliend,jend+length*normaljend,iend-length*normaliend,jend-length*normaljend],width=2,fill="grey")

        self.display.create_line([istart + il * normalistart, jstart + il * normaljstart, istart - ol * normalistart, jstart - ol * normaljstart],
                                 width=2, fill="grey")
        self.display.create_line([iend + il * normaliend, jend + il * normaljend, iend - ol * normaliend, jend - ol * normaljend],
                                 width=2, fill="grey")

    def normal_vector(self,phi):
        """Calculates a normal vector for the current segmentation for a given angle

        Args:
            phi (int): angle at which the the normal vector should be calculated

        Returns:
            normal_i (float): x-coordinate endpoint of normal vector
            normal_j (float): y-coordinate endpoint of normal vector

        """
        ileft  = fitfunc(self.cell.fit_x[self.cell.current_layer][1],phi-0.001)
        jleft  = fitfunc(self.cell.fit_y[self.cell.current_layer][1],phi-0.001)
        iright = fitfunc(self.cell.fit_x[self.cell.current_layer][1],phi+0.001)
        jright = fitfunc(self.cell.fit_y[self.cell.current_layer][1],phi+0.001)

        normal_j = -1*(ileft-iright)
        normal_i = jleft-jright
        length_normal = math.sqrt(normal_i**2+normal_j**2)
        normal_j /= length_normal
        normal_i /= length_normal

        return normal_i,normal_j
    
    def linescan(self):
        """Performs linescan analysis on selected slices. Automatically draws
        and displays linescan plots and saves the linescans to the image directory
        as .dat text files.

        Note: if pixels from the linescan fall out of the image, they are given an
        intensity value of zero.

        """
        self.store_linescan_entries()

        #picks image slices to work on
        scanlayers =  self.generate_framelist()

        for layer in scanlayers:
                  
            self.cell.seek(layer)
            self.load_linescan_entries()
            self.draw()
            self.display.update_idletasks()

            #gets linescan parameters from GUI and gets the necessary PHIs
            if self.cell.fit_x[self.cell.current_layer][0]==True:
                
                fit_x = self.cell.fit_x[self.cell.current_layer][1]
                fit_y = self.cell.fit_y[self.cell.current_layer][1]

                outer = float(self.outerlength.get())
                inner = float(self.innerlength.get())

                phistart = float(self.phistart.get())
                phiend = float(self.phiend.get())

                self.cell.outline_pixel_list[self.cell.current_layer] = []

                icurrent = int(round(fitfunc(fit_x,phiend*np.pi/180.0)))
                jcurrent = int(round(fitfunc(fit_y,phiend*np.pi/180.0)))

                self.cell.outline_pixel_list[self.cell.current_layer].append([icurrent,jcurrent,phiend*np.pi/180.0])

                for PHI in np.arange(phiend*np.pi/180.0,phistart*np.pi/180.0,-0.001):
                           
                    i = int(round(fitfunc(fit_x,PHI)))
                    j = int(round(fitfunc(fit_y,PHI)))

                    if i!=icurrent or j!=jcurrent:
                        self.cell.outline_pixel_list[self.cell.current_layer].append([icurrent,jcurrent,PHI])
                        icurrent = i
                        jcurrent = j
               
                self.draw()
                self.display.update_idletasks()

                unrolled_image = Image.new("RGB",(len(self.cell.outline_pixel_list[self.cell.current_layer]),int(outer+inner+1)),(0,0,0))
                unrolled_matrix = np.zeros((len(self.cell.outline_pixel_list[self.cell.current_layer]),int(outer+inner+1)))

                averages = np.zeros(int(inner)+int(outer)+1)
                avpixel = np.zeros((1,1))

                stepsize = 0.1
                steps_per_pixel = int(1./stepsize)  #must be divisible by 2

                outersteps = int(outer/stepsize)
                innersteps = int(inner/stepsize)
                pixels = np.zeros((len(self.cell.outline_pixel_list[self.cell.current_layer]),innersteps+outersteps+1))

                for i in range(0,len(self.cell.outline_pixel_list[self.cell.current_layer]),1):
                     
                    PHI = self.cell.outline_pixel_list[self.cell.current_layer][i][2]
                    outlinex = fitfunc(fit_x,PHI)
                    outliney = fitfunc(fit_y,PHI)
                     
                    #defines normal at the current point in the contour
                    ileft = fitfunc(fit_x,PHI+0.001)
                    jleft = fitfunc(fit_y,PHI+0.001)
                    iright = fitfunc(fit_x,PHI-0.001)
                    jright = fitfunc(fit_y,PHI-0.001)

                    normal_j = -1*(ileft-iright)
                    normal_i = jleft-jright
                    length_normal = math.sqrt(normal_i**2+normal_j**2)
                    normal_j /= length_normal
                    normal_i /= length_normal

                    current_normal_i = outlinex
                    current_normal_j = outliney
               
                    if avpixel.shape[0]==1:
                        avpixel = np.zeros(innersteps+outersteps+1)
                        counter = 0

                    center = int(len(pixels)/2)

                    #gets intensity values for inner linescan (by linear interpolation)
                    for step in range(0,innersteps,1):
                        
                        try:
                            realx = outlinex-(inner-step*stepsize)*normal_i
                            realy = outliney-(inner-step*stepsize)*normal_j

                            Q11   = int(np.ceil(realx)) ,int(np.ceil(realy))
                            Q12   = int(np.floor(realx)),int(np.ceil(realy))
                            Q22   = int(np.floor(realx)),int(np.floor(realy))
                            Q21   = int(np.ceil(realx)) ,int(np.floor(realy))

                            AQ11  = abs(realx-np.ceil(realx)) *abs(realy-np.ceil(realy))
                            AQ12  = abs(realx-np.floor(realx))*abs(realy-np.ceil(realy))
                            AQ22  = abs(realx-np.floor(realx))*abs(realy-np.floor(realy))
                            AQ21  = abs(realx-np.ceil(realx)) *abs(realy-np.floor(realy))

                            FxQ11 = np.sum(self.cell.cellimagematrix[Q11[0],Q11[1]])
                            FxQ12 = np.sum(self.cell.cellimagematrix[Q12[0],Q12[1]])
                            FxQ22 = np.sum(self.cell.cellimagematrix[Q22[0],Q22[1]])
                            FxQ21 = np.sum(self.cell.cellimagematrix[Q21[0],Q21[1]])

                            grey_level = FxQ11*AQ22+FxQ12*AQ21+FxQ22*AQ11+FxQ21*AQ12
                            pixels[i][step] = grey_level
                           
                        except IndexError:
                            print("Linescan pixel out of image")

                    #gets intensity values for inner linescan (by linear interpolation)
                    for step in range(0,outersteps+1,1):
                         
                        try:
                           realx = outlinex+(step*stepsize)*normal_i
                           realy = outliney+(step*stepsize)*normal_j

                           Q11   = int(np.ceil(realx)) ,int(np.ceil(realy))
                           Q12   = int(np.floor(realx)),int(np.ceil(realy))
                           Q22   = int(np.floor(realx)),int(np.floor(realy))
                           Q21   = int(np.ceil(realx)) ,int(np.floor(realy))

                           AQ11  = abs(realx-np.ceil(realx)) *abs(realy-np.ceil(realy))
                           AQ12  = abs(realx-np.floor(realx))*abs(realy-np.ceil(realy))
                           AQ22  = abs(realx-np.floor(realx))*abs(realy-np.floor(realy))
                           AQ21  = abs(realx-np.ceil(realx)) *abs(realy-np.floor(realy))

                           FxQ11 = np.sum(self.cell.cellimagematrix[Q11[0],Q11[1]])
                           FxQ12 = np.sum(self.cell.cellimagematrix[Q12[0],Q12[1]])
                           FxQ22 = np.sum(self.cell.cellimagematrix[Q22[0],Q22[1]])
                           FxQ21 = np.sum(self.cell.cellimagematrix[Q21[0],Q21[1]])

                           grey_level = FxQ11*AQ22+FxQ12*AQ21+FxQ22*AQ11+FxQ21*AQ12
                           pixels[i][step+innersteps] = grey_level

                        except IndexError:
                           print("Linescan pixel out of image.")

                    #fills in average linescan, unrolled image/matrix for the first index
                    av =  float(np.sum(pixels[i][0:int(steps_per_pixel/2)]))/float(steps_per_pixel/2)
                    unrolled_image.putpixel((i,0),(int(av),int(av),int(av)))
                    unrolled_matrix[i,0] = av
                    averages[0] += av

                    #fills in average linescan, unrolled image/matrix for the remaining indices
                    for block in np.arange(1,inner+outer,1):
                              
                        av = np.sum(pixels[i][int(block*steps_per_pixel-steps_per_pixel/2):int(block*steps_per_pixel+steps_per_pixel/2)])/steps_per_pixel
                        unrolled_image.putpixel((i,int(block)),(int(av),int(av),int(av)))
                        unrolled_matrix[i,int(block)]=av
                        averages[int(block)] += av
                           
                    av = float(np.sum(pixels[i][int(-steps_per_pixel/2):]))/float(steps_per_pixel/2)
                    averages[int(inner)+int(outer)] += av

                    avpixel += pixels[i]
                    counter += 1.0

                    #displays progress during linescan acquisition
                    if np.fmod(i,20)==0:
                        for j in np.arange(0,outer,0.1):
                            self.display.create_rectangle([np.floor(outlinex+j*normal_i),np.floor(outliney+j*normal_j),np.floor(outlinex+j*normal_i),np.floor(outliney+j*normal_j)],fill="blue",width=0)
                        for j in np.arange(0,inner,0.1):
                            self.display.create_rectangle([np.floor(outlinex-j*normal_i),np.floor(outliney-j*normal_j),np.floor(outlinex-j*normal_i),np.floor(outliney-j*normal_j)],fill="yellow",width=0)
                        self.display.update_idletasks()
               
            if int(self.linescanimage_toggle.get())==1:
                unrolled_image.show()
               
            #generates an average linescan
            averages /= float(len(self.cell.outline_pixel_list[self.cell.current_layer]))

            #saves average linescan
            ofile_name = self.cell.directory+"/"+self.cell.filename.split("/")[-1][:-4]+"_frame_%i_linescan_%i_%i_%i_%i_average.dat"%(self.cell.current_layer,phistart,phiend,outer,inner)
            print("Data written to %s!"%ofile_name)
            ofile = open(ofile_name,"w")
            for pixel in range(0,int(inner+outer),1):
                ofile.write("%i %f\n"%(pixel,averages[pixel]))
            ofile.close()

            #saves all linescans together
            ofile_name = self.cell.directory+"/"+self.cell.filename.split("/")[-1][:-4]+"_frame_%i_linescan_%i_%i_%i_%i_all_linescans.dat"%(self.cell.current_layer,phistart,phiend,outer,inner)
            print("Data written to %s!"%ofile_name)
            ofile = open(ofile_name,"w")
            #unrolled_matrix = unrolled_image.load()
            for pixel in range(0,int(inner+outer),1):
                ofile.write("%i "%pixel)
                for i in range(0,len(self.cell.outline_pixel_list[self.cell.current_layer]),1):
                    ofile.write("%f "%(unrolled_matrix[i,pixel]))
                ofile.write("\n")
            ofile.close()

            # #saves average linescan (high spatial resolution)
            # ofile_name = self.cell.directory+"/"+self.cell.filename.split("/")[-1][:-4]+"_frame_%i_linescan_%i_%i_%i_%i_average_fine.dat"%(self.cell.current_layer,phistart,phiend,outer,inner)
            # ofile = open(ofile_name,"w")
            # for pixel in range(0,len(avpixel),1):
            #     ofile.write("%f %f\n"%(pixel*stepsize,avpixel[pixel]/counter))
            # ofile.close()

            # #saves all linescans together (high spatial resolution)
            # ofile_name = self.cell.directory+"/"+self.cell.filename.split("/")[-1][:-4]+"_frame_%i_linescan_%i_%i_%i_%i_all_linescans_fine.dat"%(self.cell.current_layer,phistart,phiend,outer,inner)
            # print("Data written to %s!"%ofile_name)
            # ofile = open(ofile_name,"w")
            # for pixel in range(0,int(innersteps+outersteps),1):
            #     ofile.write("%f "%pixel)
            #     for i in range(0,len(self.cell.outline_pixel_list[self.cell.current_layer]),1):
            #         ofile.write("%f "%(pixels[i,pixel]))
            #     ofile.write("\n")
            # ofile.close()
               
            #displays the average linescans from the stack together
            plt.plot(range(int(-inner),int(outer)+1,1),averages)

        plt.show()

def main():
    """Starts the program by opening an application window with a default start image"""

    # master = Tk() #moved this up to the imports
    master.resizable(width=0, height=0)
    master.title(string=sys.argv[0][:-3].split("/")[-1])
    app = App(master)
    mainloop()

    #closes the program when the window is closed
    sys.exit()

if __name__=="__main__":

    main()