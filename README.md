# Cortex Thickness Analysis

Cortex thickness analysis is performed in two steps:
1) Segmentation of cortex and membrane images with `cortex_segmentor_v1_1.py`
2) Analysis of linescans with `extract_cortex_thickness_v5.py`

##  Dependencies

* [Python](www.python.org) written and tested in 2.7, updated to work on 3.X
* [Numpy](www.numpy.org)
* [Scipy](www.scipy.org)
* [Matplotlib](www.matplotlib.org)
* [PIL](http://www.pythonware.com/products/pil/)
* A window system (e.g. XQuartz for OSX)

## Additional notes

To ensure accuracy of cortex thickness measurements, your images should be color corrected both for chromatic shift and
magnification differences between colors.

This software has been used in the following publications for cell segmentation for cortex thickness/density measurements:
- **Clark A.G.**, Dierkes K. and E.K. Paluch (2013) Monitoring Actin Cortex Thickness in Live Cells. *Biophy. J.* 105(3):570-580.
- **Chugh P., Clark A.G., Smith M.B.**, Cassani D.A.D., Ragab A., Roux P.P., Charras G., Salbreux G. and E.K. Paluch (*in press*) Actin cortex architecture regulates cell surface tension. *Nat. Cell Biol.*

# Cortex Segmentor

A GUI-based program to read images, perform segmentations of an enriched region
in the cell periphery (the cortex) and generate linescans normal to the segmentation contour.

## Loading/Saving an image / segmentation

- Image stacks (tif format) can be opened using the `Open file` button
- Segmentation files (.cmd format) can be opened using the `Load` button
- Segmentation files (.cmd format) can be saved using the `Save` button


## Navigating your image stack

- You can navagate through your stack using the `Previous layer` and `Next layer` buttons.
- Addtionally, you can use the left and right arrows on your keyboard.
- You can choose what is displayed by toggling the boxes on the left.

## Cortex Segmentation

This software segments cells by first selecting pixels on the image to be considered for segmentation.
These "fit points" are then fit to a Fourier expansion to yield a smooth segmentation contour.

### Automated fit point selection

Automated fit point selection is performed by simple intensity thresholding and sequential removal/refinement of fit points.

- Fit points can be determined in an automated fashion using the `Choose fit points` button.
- Fit points will be chosen for the slices according to analysis `Mode.`
The analysis mode can be adjusted using the dropdown menu to the right of `Mode` (default:`Current`).
If the mode is `Entry`, write in a frame number in the box to work on a particular frame.
- The slider next to `Choose fit points` will determine the initial threshold used for fit point selection.
- You can watch the progress of fit point refinement if the `Fit point` and `Segmentation` boxes are ticked.
- If the automated fit point selection does not faithfully segment your data you can:
  1. choose different parameters for refining fit points in the `choose_fit_points` function of the `App` class.
  2. refine the automated segmentation by hand by manually by adding/deleting fit points
  or simply adding all of the fitpoints by hand (see next section).
- Currently, the fit point selection is optimized for segmentation of the membrane channel, with low cytoplasmic background.

### Manual fit point selection

- Fit points can be manually added and removed by right and left mouse clicks, respectively.
- The size of the fit points to be added/removed is set by the Fitpoint radius (near the bottom).
- All current fit points can be reset using the `Reset fitpoints` button.

### Segmentation

- Segmentation is accomplished by fitting the fit points with a Fourier expansion to achieve a smooth contour.
- The number of modes of the expansion can be adjusted in the text box next to the `Segment` button.
A higher number of modes will make the segmentation allows for capturing higher frequency features in the contour.
- Segmentation can be performed on a single slice or multiple slices by choosing the `Mode` (as for choosing fit points above).

### Copying segmentation data

Fit point and segmentation data can be copied from one frame to another.

- Copying is performed using the `Copy segmentation` button. This will copy fit points, the segmentation and linescan angles and lengths (see below).
- Copying can be performed in different ways using the dropdown menu next to the `Copy segmentation` button.
  - `Entry` will copy data from the frame number specified in the box next to the dropdown menu to the current frame.
  - Data can be copied from all even to odd frames and vice-versa (or by multiples of three).
  This is useful if you have images of the same cell in two (or three) channels to ensure that both channels have the exact same segmentation.

## Generating linescans

Linescans can be generated in a direction normal to the segmentation contour.

- Generation of linescans can be performed using the `Linescan` button. As for choosing fit points and segmentation,
the frames to be processed is determined by the `Mode` (same as for choosing fit points and segmenting).
- The region of the contour can be adjusted by changing the `Startangle` and `Endangle.`
To view where the start and end angles are, click the `Show endpoints` button.
- The length of the linescans can be adjusted by changing the `Outer length` and `Inner length.`
To view the linescan lengths on your image, again, click the `Show endpoints` button.
- The linescan parameters can be copied with the button `Copy to all`, which copies the current linescan parameters
(angles and lengths) to all other frames.
- When you click `Linescan,` and `Linescan - Points` is checked, you will see the progress of linescan generation on the image.
- After linescans are generated for a slice, the linescans will be saved as unformatted `.dat` files. For each image, two
linescan files are generated:
  1. An `all_linescans` file contains all of the linescans for the image, indexed by the linescan position (inner->outer position).
  2. An `average` linescan file containing the average intensity of all linescans, indexed by the linescan position (inner->outer position).
- The linescan file names contain the frame number, start/end angles and inner/outer linescan lengths.
- After the generation of linescans is finished, a plot of all of the linescans will open automatically.
- If `Linescan image` is clicked at the time of linescan generation, an image of the unrolled linescan (corresponding to `all_linescans`),
will be automatically opened.

# Cortex Thickness Extraction

The `extract_cortex_thickness` script is for the extraction of cortex thickness and density and linescan properties from linescans generated
with `cortex_segmentor`. The extraction is based on fitting a simple model of cortex geometry to paired linescans of
actin (or some other cortical marker that colocalizes with cortical actin) and the plasma membrane.

## Linescans to Analyze

The script is designed to analyze average linescans generated with `cortex_segmentor`.
These `.dat` files should contain a frame number in their name (for batch processing).

Linescans are analyzed in pairs (one for the cortical material--typically actin--and the plasma membrane).

Linescans can be analyzed as a single pair or in batch for many directories of linescan pairs.

## Extracting thickness and linescan parameters

Upon executing `extract_cortex_thickness_v5.py`, you will be prompted to continue in *single pair mode* or *batch mode*

### Single pair mode

When running the script in single pair mode, you will be further prompted to enter the following information:

1) The file path of the first linescan (channel 1)
2) The file path of the second linescan (channel 2)
3) The pixel size for your linescan (same as the original images from which linescans were generated)
4) The actin channel (channel 1 or 2)
5) The sigma value of the point spread function (PSF) of the actin channel*

*The PSF can be determined by imaging a sub-resolution object (e.g. fluorescent beads) using the same imaging settings
 as for your cortex images. The resulting signal can be fit with a 2D Gaussian to extract sigma. For details, see:
 **Clark A.G.**, Dierkes K. and E.K. Paluch (2013) Monitoring Actin Cortex Thickness in Live Cells. *Biophy. J.* 105(3):570-580.

#### Output

In single pair mode, a directory corresponding to the name of the first linescan plus `_ls_data_v5` will be made.
Plots of the linescans with and without fits of the cortex model and a data list with the resulting cortex measurements
and linescan properties called `ls_data.dat`.

### Batch mode

When running the script in batch mode, you will be prompted to choose a `parent_directory` that contains sub-directories
with linescans and a `dir_list.dat` file containing a list of the sub-directories to be analyzed and some required parameters.
As noted in the docstrings, `dir_list.dat` file should have the following format:

```
sub_dir  px_size  category  ch_actin  sigma_actin
stk_1    0.05     control   1         0.119
stk_2    0.04     siRNA     2         0.220
...
```

With the following definitions of the parameters:
- `sub_dir`: The name of the sub-directory containing the linescan pairs
- `px_size`: The pixel size for the linescans in the given sub_dir
- `category`: The category of the experiment in each sub_dir (used for keeping track of different experiments)
- `ch_actin`: The actin channel (either 1 or 2)
- `sigma_actin`: The sigma of the point spread function for the actin channel (as described above)

Linescans will be paired according to frame numbers, which should be given in the linescan filenames
(e.g., `frame_0`/`frame_1`, `frame_2`/`frame_3`, etc).

#### Output

In batch mode, in each sub-directory, a directory called `_ls_data_v5` will be made.
Plots of the linescans with and without fits of the cortex model and a data list with the resulting cortex measurements
and linescan properties called `ls_data.dat`. The data rows are labeled according to the name of the channel 1 linescan.

In addition, a text file called `master_list_v5.dat` will be created in the `parent_dir` that combines all of the analyzed
data together.

# Known issues and disclaimers

- The `cortex_segmentor` GUI was designed to handle 512x512 pixel images. You may run into problems if you images are not this size.
If your images are smaller, you can make them 512x512 by adding a black border around using e.g. FIJI.
- The fits for the segmentation implies a center near the center of the image. You will have trouble segmenting your
cells if this is not the case.

# TODO

- Allow for window resizing and handling different sized images for `cortex_segmentor`
- Rearrange layout of GUI to improve organization
- Allow for easier modification of segmentation parameters (and saving/loading segmentation parameters) for `cortex_segmentor`
- Improve mode/parameter selection via user interface for `extract_cortex_thickness`
- Combine `extract cortex_thickness` with `cortex_segmentor` to integrate segmentation and thickness extraction
