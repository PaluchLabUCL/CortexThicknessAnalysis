#!/opt/local/bin/python

__author__ = "Andrew G. Clark"
__date__ = "7 May 2014"
__maintainer__ = "Andrew G. Clark"
__email__ = "andrew.clark@curie.fr"

""" This script analyzes linescans and extracts cortex thickness and density from actin/membrane linescan pairs.

The script can be run in a 'pair' mode (to analyze a single linescan pair)
or 'batch' mode (to analyze multiple directories full of linescan pairs).
The mode can be specified at the bottom ("main" function).

For batch mode:

Your parent directory should contain a file called 'dir_list.dat'
with the following information in row/column form, with only space as delimiters:

sub_dir  px_size  category  ch_actin  sigma_actin
stk_1    0.05     control   1         0.119
stk_2    0.04     siRNA     2         0.220
...

The first row must contain the column headers as shown
Definitions of input parameters:

sub_dir: The name of the sub-directory containing the linescan pairs (linescan pairs must end in '...average.dat')
px_size: The pixel size for the linescans in the given sub_dir
category: The category of the experiment in each sub_dir (can be used for plotting later)
ch_actin: The actin channel (either '1' or '2'; used for extracting cortex thickness/i_c)
sigma_actin: The sigma of the point spread function for the actin channel (used for extracting h/i_c)

Note: For the sub_dir entries in the dir_list, only those directories NOT appearing in 'completed_list_v4_1.dat' will be analyzed

Output:

In each sub-directory, a list called '.../ls_data/ls_fit_data.dat' will be created containing linescan and thickness data
    -The columns are labeled according to channel number (ch1/ch2)
    -delta is always the position of the peak intensity of channel 2 (ch2.x_peak) minus ch1.x_peak
In each sub-directory, plots of the linescans and the linescans with fits (if applicable) will be saved in '.../ls_plots/'
At the end, a master list of all of the data combined is be created in the parent_directory

For 'manual' mode:

When running the script, windows will pop up sequentially to request the following information:

-Channel 1 average linescan file
-Channel 2 average linescan file
-Pixel Size
-Actin Channel
-Sigma (Actin)

These parameters are defined above.

"""

import os
import math
from copy import deepcopy
from tkinter import *
from tkinter.filedialog import *
from tkinter.simpledialog import *

root = Tk()

import scipy
from scipy import optimize, stats
import pylab
import numpy as np

import utility_functions as uf

def gauss_func(p, x):
    """Definition of gaussian function used to fit linescan peaks.
    p = [a, sigma, mu, c].

    """
    a, sigma, mu, c = p #unpacks p (for readability)
    g = a / (sigma * math.sqrt(2 * math.pi)) * scipy.exp(-(x - mu)**2 / (2 * sigma**2)) + c
    return g

def convolved(p,x):
    """Defines convolved linescan. Args: x: float or list/iterable of floats,
    the position for which convolved intensity is calculated; p: list/iterable
    of floats, linecan parameters (p=[i_in, i_c, i_out, h, x_c, sigma]).
    Returns: i: float, intensity at x.

    """
    i_in, i_c, i_out, h, x_c, sigma = p #unpacks p (for readability)

    i = (i_in + (i_c - i_in) * stats.norm.cdf((x - x_c) + h / 2., 0., sigma) +
         (i_out - i_c) * stats.norm.cdf((x - x_c) - h / 2., 0., sigma))

    return i

def unconvolved(p,x):
    """Defines unconvolved linescan. Args: x: float or list/iterable of floats,
    the position for which intensity is calculated; p: list/iterable of floats,
    linecan parameters (p=[i_in, i_c, i_out, h, x_c]). Returns: i: float,
    intensity at x.

    """

    i_in, i_c, i_out, h, x_c = p #unpacks p (for readability)

    i = np.zeros(len(x))

    for j in range(len(x)):
        if x[j] < x_c - h / 2.:
            i[j] = i_in
        if x[j] >=  x_c - h / 2. and x[j] <  x_c + h / 2.:
            i[j] = i_c
        if x[j] >= x_c + h / 2.:
            i[j] = i_out

    return i

def sort_ls_list(list):
    """Sorts list of linescan files by keyword.

    Args:
        list (list): the list to be sorted (here, linescan filenames)
        param (str): the keyword to use for sorting (here, usually 'frame')

    """

    def find_key(line):
        key = int(re.search('frame_([0-9]+)_', line).group(1))
        return key

    list.sort(key=find_key)
    return list

class Linescan():
    """Linescan object with methods to extract important parameters
    from linescans.

    """

    def __init__(self,x,i):
        """Initializes linescan.

        Args:
            x (list of numbers): the position values
            i (list of numbers): the intensity values

        """
        #populate linescan position/intensity
        self.x = np.array(x,dtype='float') #position list as NumPy array of floats
        self.i = np.array(i,dtype='float') #intensity list as NumPy array of floats

        #detminere a few easy parameters from position/intensity
        self.H = self.x[-1] - self.x[0]
        self.i_tot = np.trapz(self.i,self.x)

        #populate other attributes
        self.dist_to_x_in_out = 1. #specifies how far away x_in is from the peak (in um)
        self.gauss_params = None #parameter list from Gaussian fit to find peak
        self.x_peak = None #linescan peak position
        self.i_peak = None #linescan peak intensity
        self.i_in = None #intracellular intensity
        self.i_out = None #extracellular intensity
        self.max_idx = None #index of point near linescan center with highest intensity
        self.x_fit = None #position list used for peak fitting
        self.i_fit = None #intensity list used for peak fitting
        self.i_in_x_list = None #position list used to determine self.i_in
        self.i_in_i_list = None #intensity list used to determine self.i_in
        self.i_out_x_list = None #position list used to determine self.i_out
        self.i_out_i_list = None #intensity list used to determine self.i_out
        self.x_in_upper_index = None #the index at the upper end of the region where x_in is calculated
        self.x_out_lower_index = None #the index at the lower end of the region where x_out is calculated
        self.fwhm = None #full width at half-max

        #initializes linescans and determines linescan parameters
        self.extract_ls_parameters()

    def convert_px_to_um(self):
        """Multiplies list of coordinates by pixel_size."""

        self.x = np.array([a * self.px_size for a in self.x])

    def extract_ls_parameters(self):
        """Extracts intensity and position information from linescan"""

        self.get_peak()
        self.get_i_in_out()
        self.get_fwhm()

    def get_peak(self):
        """Finds the peak position and intensity of a linescan by fitting
        a Gaussian near the peak.

        """

        #restricts fitting to near the center of the linescan
        self.max_idx = int(np.argmax(self.i[int(len(self.i)/2-6):int(len(self.i)/2+20)]) + len(self.i)/2-6)
        self.x_fit = self.x[int(self.max_idx-2):int(self.max_idx+3)]
        self.i_fit = self.i[int(self.max_idx-2):int(self.max_idx+3)]

        #picks reasonable starting values for fit
        self.i_in_guess = np.mean(self.i[:int(self.max_idx-14)])
        a = (self.i[self.max_idx] - self.i_in_guess) / 2.4
        sigma = 0.170
        mu = self.x[self.max_idx]
        b = self.i_in_guess

        #perform fit with starting values
        p0 = [a, sigma, mu, b]
        p1, success  = optimize.leastsq(self.residuals_gauss,p0,
                                        args=(self.x_fit, self.i_fit),
                                        maxfev = 1000000)
        self.gauss_params = p1
        self.x_peak = p1[2]
        self.i_peak = gauss_func(p1, self.x_peak)

    def get_i_in_out(self):
        """Gets values for intracellular intensity (self.i_in) and
        extracellular intensity (self.i_out). The left of the linescan
        (nearer zero) is always assumed to be the intracellular side.
        Note: the i_in and i_out values are calculated to be the average value
        of the ten points out from the distance between the peak and position x away
        from the peak, where x is given by self.dist_to_x_in_out (defined in __init__).
        """

        x_in_upper = self.x_peak - self.dist_to_x_in_out
        x_in_upper_index = np.argmin(abs(self.x - x_in_upper))
        self.x_in_upper_index = x_in_upper_index #for use in finding total intensity for density calculation
        self.i_in_x_list = self.x[x_in_upper_index-10:x_in_upper_index]
        self.i_in_i_list = self.i[x_in_upper_index-10:x_in_upper_index]
        self.i_in = np.mean(self.i_in_i_list)

        x_out_lower = self.x_peak + self.dist_to_x_in_out
        x_out_lower_index = np.argmin(abs(self.x - x_out_lower))
        self.x_out_lower_index = x_out_lower_index #for use in finding total intensity for density calculation
        self.i_out_x_list = self.x[x_out_lower_index:x_out_lower_index+10]
        self.i_out_i_list = self.i[x_out_lower_index:x_out_lower_index+10]
        self.i_out = np.mean(self.i_out_i_list)

    def residuals_gauss(self,p,x,x_data):
        """Returns residuals for Gaussian fit of the intensity peak.
        Possible values for fit parameters are constrained to avoid
        overestimation of peak intensity.

        Args:
            p (list): fit parameters, [a, sigma, mu, c]
            x (list): position values
            x_data (list): intensity values

        Returns:
            residuals (list): residuals for fit
             -or-
            fail_array (list): in place of residuals if the fit fails

        """

        a, sigma, mu, c = p #unpacks p (for readability)

        i_peak_guess = gauss_func(p, mu)

        fail_array = np.ones(len(x)) * 99999.

        if all([sigma >= 0.1,
               abs(i_peak_guess - self.i[self.max_idx]) < 0.5 * self.i[self.max_idx]]):

            residuals = gauss_func(p,x) - x_data
            return residuals

        else:
            return fail_array

    def get_fwhm(self):
        """Calculates the full-width at half maximum (FWHM) of the linescan peak"""

        #determines half-max
        hm = (self.i_in + self.i_peak) / 2.
        # print(hm)

        # finds points closest to hm to the left of the peak
        search = self.i[:self.max_idx]
        self.left_index = (np.abs(search - hm)).argmin()
        if hm > self.i[self.left_index]:
            self.left_index_left = deepcopy(self.left_index)
            self.left_index_right = self.left_index_left + 1
        else:
            self.left_index_right = deepcopy(self.left_index)
            self.left_index_left = self.left_index_right - 1

        #gets interpolated intensity (linear interpolation between 2 surrounding points
        m_left = (self.i[self.left_index_right] - self.i[self.left_index_left]) /  (self.x[self.left_index_right] - self.x[self.left_index_left])
        b_left = self.i[self.left_index_right] - m_left * self.x[self.left_index_right]
        x_fwhm_left = (hm - b_left) / m_left
        self.fwhm_left = [x_fwhm_left,hm]

        #finds point closest to hm to the right of the peak
        search = self.i[self.max_idx:]
        self.right_index = (np.abs(search - hm)).argmin() + self.max_idx
        if hm < self.i[self.right_index]:
            self.right_index_left = deepcopy(self.right_index)
            self.right_index_right = self.right_index_left + 1
        else:
            self.right_index_right = deepcopy(self.right_index)
            self.right_index_left = self.right_index_right - 1

        #gets interpolated intensity (linear interpolation between 2 surrounding points
        m_right = (self.i[self.right_index_right] - self.i[self.right_index_left]) / (self.x[self.right_index_right] - self.x[self.right_index_left])
        b_right = self.i[self.right_index_right] - m_right * self.x[self.right_index_right]
        x_fwhm_right = (hm - b_right) / m_right
        self.fwhm_right = [x_fwhm_right,hm]

        self.fwhm = x_fwhm_right - x_fwhm_left

class Cortex():
    """A Class for a cortex, with actin and membrane linescans and
     methods to determine cortex thickness and density.


    """
    def __init__(self,ch1,ch2,sigma_actin,ch_actin=1):
        """Initializes linescan pairs and remaining attributes.

            Args:
                ch1 (Linescan class): the ch1 linescan
                ch2 (Linescan class): the ch2 linescan
                sigma_actin (float): the sigma of the PSF for the actin channel

            Kwargs:
                ch_actin (int): says which channel is actin

        """
        self.ch1 = ch1
        self.ch2 = ch2
        self.sigma_actin = sigma_actin
        self.ch_actin = ch_actin

        self.delta = self.ch2.x_peak - self.ch1.x_peak #separation between ch2 and ch1 peaks

        if self.ch_actin==1:
            self.actin = self.ch1
            self.memb = self.ch2
        elif self.ch_actin==2:
            self.actin = self.ch2
            self.memb = self.ch1
        else:
            self.actin = None
            self.memb = None

        self.h_max = 1. #maximum cortex thickness (for constraining fit)
        self.i_c_max = 500. #maximum cortex intensity (for constraining fit)
        self.h = None #cortex thickness (from fit)
        self.i_c = None #cortical actin intensity (from fit)
        self.density = None #cortical actin density
        self.X_c = None #background-independent center position of the cortical actin (from fit)
        self.solution = None #solution from actin cortex thickness fit

    def get_h_i_c(self):
        """ Performs fit to get cortex thickness, h, and cortex intensity, i_c

         Note: density is calculated as the difference between fitted cortex intensity
         and intracellular background, normalized by the intensity from the beginning
         of the linescan to end of the i_out calculation region

        """

        delta = abs(self.delta)

        #SET STARTING VALUES FOR ROOTS AND SOLUTIONS
        self.solution = 2e+20

        #only try fitting if the peak is higher than both i_in and i_out
        if ((self.actin.i_out - self.actin.i_peak) /
                (self.actin.i_in - self.actin.i_peak))>=0:

            #loops through several different starting values for i_c and h
            for i_c_factor in np.arange(2.,3.1,0.2):
                for h_factor in np.arange(0.5, 2.1, 0.2):

                    i_c_start = self.actin.i_peak * i_c_factor
                    delta_start = ((self.sigma_actin**2 / delta*2) *
                                   np.log((self.actin.i_out - i_c_start) /
                                          (self.actin.i_in - i_c_start)))
                    h_start = 2 * (delta - delta_start) * h_factor

                    #performs fit
                    p0 = [h_start, i_c_start]

                    try:
                        result = optimize.leastsq(self.residuals, p0,
                                                  maxfev=100000, full_output=1)

                        solution_temp = np.sum([x**2 for x in result[2]['fvec']])

                        if solution_temp < self.solution:
                            self.solution = deepcopy(solution_temp)
                            p1 = result[0]

                    except TypeError:
                        pass

            #controls for bad fits
            if any([self.solution>0.01,
                    p1[0] >= self.h_max - 0.001,
                    p1[1] >= self.i_c_max - 1.]):
                 p1 = [None, None]
                 self.h = None
                 self.i_c = None
                 self.density = None
                 self.X_c = None
                 self.solution = None
            else:
                self.h, self.i_c = p1
                actin_ls_mean = np.mean(self.actin.i[:self.actin.x_out_lower_index+10])
                self.density = (self.i_c - self.actin.i_in) / actin_ls_mean
                self.X_c = self.memb.x_peak - self.h / 2.

    def residuals(self,p):
        """Calculates residuals for cortex linescan fit to extract cortex
        thickness and intensity values

        Args:
            p (list of floats): [thickness, cortex_intensity]

        Returns:
            residuals (list of floats): [residual1, residual2]
            -or-
            fail_array (list of floats): [1000000., 1000000.]
             (returned only if fitting fails)

        """

        fail_array = [1000000., 1000000.]

        #constrains fit and ensures log term is positive
        if all([self.h_max>p[0]>0,
               self.i_c_max>p[1]>self.actin.i_in,
               (self.actin.i_out - p[1]) / (self.actin.i_in - p[1]) > 0]):

            #X_c is the position of the center of the cortex
            #x_c is the position of the cortex peak
            X_c_try = self.memb.x_peak - p[0] / 2.
            delta_try = (self.sigma_actin**2 / p[0]) * np.log((self.actin.i_out - p[1]) / (self.actin.i_in - p[1]))
            x_c_try = X_c_try - delta_try
            i_peak_try = convolved([self.actin.i_in, p[1], self.actin.i_out, p[0], X_c_try, self.sigma_actin], x_c_try)

            #residuals are difference between calculated peak position/intensity and values from data
            residuals = [x_c_try - self.actin.x_peak, i_peak_try - self.actin.i_peak]
            return residuals

        else:
            return fail_array

    def plot_lss(self):
        """Plots linescans"""

        fig = pylab.figure()
        ax = fig.add_subplot(1,1,1)

        #plots raw data
        pylab.plot(self.ch1.x,self.ch1.i,'go',label="Ch. 1")
        pylab.plot(self.ch2.x,self.ch2.i,'ro',label="Ch. 2")

        #plots points used for determining i_in and i_out
        pylab.plot(self.ch1.i_in_x_list,self.ch1.i_in_i_list,'yo',label=r"$i_{\rm{in}}$, $i_{\rm{out}}$")
        pylab.plot(self.ch2.i_in_x_list,self.ch2.i_in_i_list,'yo')
        pylab.plot(self.ch1.i_out_x_list,self.ch1.i_out_i_list,'yo')
        pylab.plot(self.ch2.i_out_x_list,self.ch2.i_out_i_list,'yo')

        #plots points used to calculate fwhm and shows the fwhm
        # pylab.plot(self.ch1.x[self.ch1.left_index_left],self.ch1.i[self.ch1.left_index_left],'ko',label="fwhm points")
        # pylab.plot(self.ch1.x[self.ch1.left_index_left],self.ch1.i[self.ch1.left_index_left],'ko')
        # pylab.plot(self.ch1.x[self.ch1.left_index_right],self.ch1.i[self.ch1.left_index_right],'ko')
        # pylab.plot(self.ch1.x[self.ch1.right_index_left],self.ch1.i[self.ch1.right_index_left],'ko')
        # pylab.plot(self.ch1.x[self.ch1.right_index_right],self.ch1.i[self.ch1.right_index_right],'ko')
        #
        # pylab.plot(self.ch2.x[self.ch2.left_index_left],self.ch2.i[self.ch2.left_index_left],'ko')
        # pylab.plot(self.ch2.x[self.ch2.left_index_right],self.ch2.i[self.ch2.left_index_right],'ko')
        # pylab.plot(self.ch2.x[self.ch2.right_index_left],self.ch2.i[self.ch2.right_index_left],'ko')
        # pylab.plot(self.ch2.x[self.ch2.right_index_right],self.ch2.i[self.ch2.right_index_right],'ko')

        x_fwhm1, i_fwhm1 = zip(self.ch1.fwhm_left,self.ch1.fwhm_right)
        x_fwhm2, i_fwhm2 = zip(self.ch2.fwhm_left,self.ch2.fwhm_right)

        pylab.plot(x_fwhm1, i_fwhm1,'g',ls='-',marker='x',label="fwhm")
        pylab.plot(x_fwhm2, i_fwhm2,'r',ls='-',marker='x',label='fwhm')

        # x_fwhm1 = [self.ch1.x[self.ch1.left_index],self.ch1.x[self.ch1.right_index]]
        # y_fwhm1 = (self.ch1.i[self.ch1.left_index] + self.ch1.i[self.ch1.right_index]) / 2.
        # i_fwhm1 = [y_fwhm1,y_fwhm1]
        # pylab.plot(x_fwhm1,i_fwhm1,'g-',label="fwhm")
        #
        # x_fwhm2 = [self.ch2.x[self.ch2.left_index],self.ch2.x[self.ch2.right_index]]
        # y_fwhm2 = (self.ch2.i[self.ch2.left_index] + self.ch2.i[self.ch2.right_index]) / 2.
        # i_fwhm2 = [y_fwhm2,y_fwhm2]
        # pylab.plot(x_fwhm2,i_fwhm2,'r-',label="fwhm")

        #plots gaussian fit curve
        x_gauss_fit_ch1 = np.linspace(self.ch1.x_fit[0],self.ch1.x_fit[-1],100)
        i_gauss_fit_ch1 = gauss_func(self.ch1.gauss_params,x_gauss_fit_ch1)
        pylab.plot(x_gauss_fit_ch1,i_gauss_fit_ch1,'b',label="Peak fit")

        x_gauss_fit_ch2 = np.linspace(self.ch2.x_fit[0],self.ch2.x_fit[-1],100)
        i_gauss_fit_ch2 = gauss_func(self.ch2.gauss_params,x_gauss_fit_ch2)
        pylab.plot(x_gauss_fit_ch2,i_gauss_fit_ch2,'b')

        #finish plot
        y_min, y_max = ax.get_ylim()
        pylab.ylim = (0,y_max)

        pylab.xlabel("Position ($\mu$m)")
        pylab.ylabel("Intensity (AU)")
        pylab.legend(loc='upper right')
        pylab.gcf().subplots_adjust(bottom=0.15)

    def plot_fits(self):
        """Plots linescan pair with fitted cortex thickness"""

        fig = pylab.figure()
        ax = fig.add_subplot(1,1,1)

        if self.ch_actin==1 or self.ch_actin=="1":
            color_actin = 'g'
            color_memb = 'r'
        elif self.ch_actin==2 or self.ch_actin=="2":
            color_actin = 'r'
            color_memb = 'g'
        else:
            raise ValueError("Please specify ch_actin as <<1>>, <<2>> for plotting fit!")

        #plots raw data
        pylab.plot(self.memb.x,self.memb.i,'o',color=color_memb,label="Memb. (raw)")
        pylab.plot(self.actin.x,self.actin.i,'o',color=color_actin,label="Actin (raw)")

        #plots unconvolved and extracted actin linescans from fits
        x_actin_hd = np.linspace(self.actin.x[0],self.actin.x[-1],1000)
        i_actin_unconv = unconvolved([self.actin.i_in, self.i_c,
                                       self.actin.i_out, self.h, self.X_c],
                                      x_actin_hd)
        i_actin_conv = convolved([self.actin.i_in, self.i_c,
                                   self.actin.i_out, self.h, self.X_c, self.sigma_actin],
                                  x_actin_hd)

        pylab.plot(x_actin_hd,i_actin_unconv,ls='-',color=color_actin, label='fit')
        pylab.plot(x_actin_hd,i_actin_conv,ls='--',color=color_actin, label='fit (conv.)')

        pylab.axvline(x=self.memb.x_peak, color=color_memb, ls='--', label="Memb. (peak)")

        #finishes plot
        y_min, y_max = ax.get_ylim()
        pylab.ylim = (0,y_max)

        pylab.xlabel("Position ($\mu$m)")
        pylab.ylabel("Intensity (AU)")
        pylab.legend(loc='upper right')
        pylab.gcf().subplots_adjust(bottom=0.15)

def write_master_list(parent_dir,version):
    """Writes a master data lis in the parent directory for batch mode.

    Args:
        parent_dir (string): path of the parent directory
        version (string): the version of the software (for naming output file)

    """

    dir_list_path = parent_dir + '/dir_list.dat'
    subdir_list = [_[0] for _ in uf.read_file(dir_list_path)][1:]

    master_data = []
    for i in range(len(subdir_list)):
        data_dir = parent_dir + '/' + subdir_list[i]
        data = uf.read_file(data_dir + '/ls_data/ls_data.dat')
        if i==0:
            for line in data:
                master_data.append(line)
        else:
            for line in data[1:]:
                master_data.append(line)

    # print master_data
    uf.save_data_array(master_data, parent_dir + '/master_list_v%s.dat'%version)

def load_ls(ls_path,px_size=1.):
    """Loads a linescan file

    Args:
        ls_path (str): path of the average linescan file to be loaded
        px_size (float): pixel size in microns

    Returns:
        x (numpy array): the positions (in microns)
        i (numpy array): the intensities

    """

    ls_data = uf.read_file(ls_path)
    x = np.array([float(_[0]) for _ in ls_data]) * px_size
    i = np.array([float(_[1]) for _ in ls_data])
    return x,i

def analyze_cortex(file_ch1,file_ch2,px_size,ch_actin,sigma_actin):

    """Extracts linescan parameters and coretx thickness/density
    for a pair of linescans

    Args:
        file_ch1 (str): the filepath for the first linescan
        file_ch2 (str): the filepath for the second linescan
        px_size (float): the pixel size for the linescans (for the whole directory)
        ch_actin (int): the channel of the actin linescan (1 or 2)
        sigma_actin (float): the sigma of the PSF for the actin channel

    Kwargs:
        category (str): used to keep track of different conditions in the output data file

    Returns:
        cortex (Cortex class): the cortex with associated attributes

    """

    x_ch1, i_ch1 = load_ls(file_ch1,px_size=px_size)
    x_ch2, i_ch2 = load_ls(file_ch2,px_size=px_size)
    x = deepcopy(x_ch1) #the x values should be the same for both linescans!

    basename = file_ch1.split('/')[-1][:-4]
    print('Analyzing file pair for:', basename)

    # extracts data
    actin = Linescan(x,i_ch1)
    memb = Linescan(x,i_ch2)
    cortex = Cortex(actin, memb, sigma_actin, ch_actin=ch_actin)

    if ch_actin==1 or ch_actin==2:
        cortex.get_h_i_c()
    elif ch_actin == "None":
        pass
    else:
        raise ValueError("Please specify ch_actin as <<1>> or <<2>> for %s!"%file_ch1)

    print('h =', cortex.h)
    return cortex

def analyze_ls_pair(file_ch1,file_ch2,px_size,ch_actin,sigma_actin,version):
    """Analyzes linescans to extract cortex thickness/density
    for a single linescan pair. Data and plots are generated and saved
    to a new folder with same name as file_ch1

    Args:
        file_ch1 (str): the filepath for the first linescan
        file_ch2 (str): the filepath for the second linescan
        px_size (float): the pixel size for the linescans (for the whole directory)
        ch_actin (int): the channel of the actin linescan (1 or 2)
        sigma_actin (float): the sigma of the PSF for the actin channel

    """

    # makes directory in data_dir for saving
    save_dir = file_ch1[:-4] + '_ls_data'
    uf.make_dir(save_dir)

    # makes a list of parameters to extract from cortex data
    data_to_write = [['basename', 'category',
                      'delta', 'h', 'i_c', 'density', 'X_c', 'solution',
                      'ch1.i_tot', 'ch1.H', 'ch1.x_peak', 'ch1.i_peak', 'ch1.i_in', 'ch1.i_out', 'ch1.fwhm',
                      'ch2.i_tot', 'ch2.H', 'ch2.x_peak', 'ch2.i_peak', 'ch2.i_in', 'ch2.i_out', 'ch2.fwhm'
                      ]]

    basename = file_ch1.split('/')[-1][:-4]
    category = 'pair'

    #gets cortex and linescan data
    cortex = analyze_cortex(file_ch1, file_ch2, px_size, ch_actin, sigma_actin)

    # plots raw linescans
    cortex.plot_lss()
    pylab.savefig(save_dir + "/" + basename + ".png")
    pylab.close()

    # plots linescans with h fits
    if cortex.h != None:
        cortex.plot_fits()
        pylab.savefig(save_dir + "/" + basename + "_fit.png")
        pylab.close()

    # gets extracted linescan data
    data_temp = [basename, category]
    for param in data_to_write[0][2:]:
        data_temp.append(eval("cortex.%s" % param))
    data_to_write.append(data_temp)

    # print data_to_write
    uf.save_data_array(data_to_write, save_dir + "/ls_data.dat")

def analyze_dir(data_dir,px_size,category,ch_actin,sigma_actin,version):
    """ Analyzes all linescan pairs in a directory full of linescans

    Args:
        data_dir (str): the directory containing the linescans
        px_size (float): the pixel size for the linescans (for the whole directory)
        category (str): the category for the experiment
        ch_actin (int): the channel of the actin linescan (1 or 2)
        version (str): version number (for output filenames)

    """

    #makes necessary directories in data_dir for saving
    save_dir = data_dir + '/ls_data'
    uf.make_dir(save_dir)

    #makes a list of parameters to extract from cortex data
    data_to_write = [['basename','category',
                      'delta', 'h', 'i_c', 'density', 'X_c', 'solution',
                      'ch1.i_tot','ch1.H','ch1.x_peak','ch1.i_peak','ch1.i_in','ch1.i_out','ch1.fwhm',
                      'ch2.i_tot','ch2.H','ch2.x_peak','ch2.i_peak','ch2.i_in','ch2.i_out','ch2.fwhm'
                      ]]

    #gets and sorts list of average linescans
    linescan_list = [x for x in os.listdir(data_dir) if 'average.dat' in x]

    for _ in linescan_list:
        print(_)
        print(re.search('frame' + '_([0-9]+)_', _).group(1))
    linescan_list = sort_ls_list(linescan_list)


    #extracts linescan parameters and thickness/density
    for i in range(int(len(linescan_list)/2)):

        file_ch1 = data_dir + '/' + linescan_list[2*i]
        file_ch2 = data_dir + '/' + linescan_list[2*i + 1]
        basename = file_ch1.split('/')[-1][:-4]

        cortex = analyze_cortex(file_ch1,file_ch2,px_size,ch_actin,sigma_actin)

        # plots raw linescans
        cortex.plot_lss()
        pylab.savefig(save_dir + "/" + basename + ".png")
        pylab.close()

        # plots linescans with h fits
        if cortex.h != None:
            cortex.plot_fits()
            pylab.savefig(save_dir + "/" + basename + "_fit.png")
            pylab.close()

        # gets extracted linescan data
        data_temp = [basename,category]
        for param in data_to_write[0][2:]:
            data_temp.append(eval("cortex.%s"%param))
        data_to_write.append(data_temp)

    # print data_to_write
    uf.save_data_array(data_to_write,save_dir + "/ls_data.dat")


def main():
    """__main__ function"""

    version = '5'

    #set up root for asking questions
    # root = Tk() #moved this up to the imports
    root.withdraw()

    #chooses analysis mode
    mode = askinteger(title="Analysis Mode Selection",
                      prompt="Please enter:\n1 for pairwise analysis or\n2 for batch analysis",
                      minvalue=1,maxvalue=2)

    if mode==1:

        ch1_path = askopenfilename(title='Select an average linescan file for channel 1',
                                   filetypes=[("dat", "*.dat")],
                                   initialdir='.',
                                   initialfile="")

        ch2_path = askopenfilename(title='Select an average linescan file for channel 2',
                                   filetypes=[("dat", "*.dat")],
                                   initialdir='/'.join(ch1_path.split('/')[:-1]),
                                   initialfile=ch1_path.split('/')[-1])

        px_size = askfloat(title='Pixel Size',prompt='Please enter your pixel size')
        ch_actin = askinteger(title='Actin Channel',prompt='Please enter the actin channel',
                              minvalue=1, maxvalue=2)
        sigma_actin = askfloat(title='Actin Sigma',prompt='Please enter the sigma value\nfor the PSF for the actin channel\n(in microns)')

        analyze_ls_pair(ch1_path,ch2_path,px_size,ch_actin,sigma_actin,version)

    if mode==2:

        parent_dir = askdirectory(title='Select the parent directory (be sure it contains dir_list.dat!)',
                                  initialdir=os.path.split(os.path.realpath(__file__))[0])
        # parent_dir = './test_data'
        dir_list = uf.get_dict_list(uf.read_file(parent_dir + '/dir_list.dat'))

        for line in dir_list:

            sub_dir = line['sub_dir']
            px_size = float(line['px_size'])
            category = line['category']
            ch_actin = int(line['ch_actin'])
            sigma_actin = float(line['sigma_actin'])
            data_dir = parent_dir + '/' + sub_dir

            print(data_dir)

            analyze_dir(data_dir,px_size,category,ch_actin,sigma_actin,version)

        write_master_list(parent_dir,version)

if __name__ == '__main__':
    main()
