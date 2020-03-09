import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.animation import FuncAnimation
from itertools import repeat, chain
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib as mpl

# random.seed( 30 )

''' import a fasta file
This is how it was created:
tail -n 100 genome.fa | head -n 200 > genome_part.txt

'''

# cdict = {'red':
# [[1.0,  1.0, 1.0],
# [1.0,  1.0, 1.0],
# [1.0,  1.0, 1.0]],
# 'blue':   [[0.0,  0.0, 0.0],
# [0.5,  1.0, 1.0],
# [1.0,  1.0, 1.0]],
# 'green': [[0.0,  0.0, 0.0],
# [0.25, 0.0, 0.0],
# [0.75, 1.0, 1.0],
# [1.0,  1.0, 1.0]],
# 'yellow':  [[0.0,  0.0, 0.0],
# [0.5,  0.0, 0.0],
# [1.0,  1.0, 1.0]],
# 'orange':  [[0.0,  0.0, 0.0],
# [0.9,  0.0, 0.0],
# [1.0,  1.0, 1.0]]}
# newcmp = LinearSegmentedColormap('testCmap', segmentdata=cdict, N=5)

cmap = plt.cm.jet  # define the colormap
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
cmaplist[0] = (0.0, 0.0, 0.0, 1.0)
cmaplist[1] = (1, 1, .5, 1.0) ## A; correct
cmaplist[2] = (.5, .5, .5, 1.0) ## C; correct
cmaplist[3] = (0.20, 1.0, 0.176, 1.0) ## G
cmaplist[4] = (1.0, 1.0, 1.0, 1.0) ## T

# create the new map
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cmap.N)

fo = open("/Users/morril01/Documents/PhD/CDA_in_Cancer/data/genome_part.txt", "r")
genome = ''.join([x.strip('\n') for x in fo.readlines()])
genome = genome.strip('n')
genome = "CTGA"*60

max_width = 500
max_height = 500
scale_timing = 40
length_read = 4*scale_timing
initial_position = random.randrange(0, len(genome)-length_read)
proposed_radius = 20
max_radius = proposed_radius

genome = "".join(chain(*zip(*repeat(genome, scale_timing))))
# max_width = 512
# max_height = 512


# Some global variables to define the whole run
total_number_of_frames = length_read
# all_data = [
#     np.random.rand(max_width, max_height) for x in range(100)
# ]



def give_int(nucl):
    if nucl == 'A':
        return (1)
    elif nucl == 'C':
        return (2)
    elif nucl == 'G':
        return (3)
    elif nucl == 'T':
        return (4)

def give_df(all_data, x, y, label_as_int):
    for j in range(y-proposed_radius, y+proposed_radius):
        if (j> 0) & (j < max_height):
            for i in range(x-proposed_radius, x+proposed_radius):
                ## if we are at a certain radius fron the centre
                if (i > 0) & (i < max_width):
                    if np.sqrt((x-i)**2 + (y-j)**2) < max_radius:
                        all_data[i,j] = label_as_int
    return all_data.tolist()

def change_val(df, label_as_int):
    return_df = np.where(df==1, label_as_int, df)
    return return_df 


def place_random_light():
    ## choose a random position in the grid
    x = random.randrange(0, max_width)
    y = random.randrange(0, max_height)
    return x,y

    ## modify all of the surrounding cells

def create_random_df():
    x,y = place_random_light()
    all_data = np.zeros(shape=(max_width, max_height))
    all_data[x,y] = 1
    # all_data = [np.array(give_df(all_data, x, y, give_int(genome[initial_position+t].upper()))) for t in range(length_read)]
    all_data = give_df(all_data, x, y, 1)
    return all_data

# def add_noise(df):
#     print( np.where(df==1))



first_df = np.array(create_random_df())

print(genome[initial_position:(initial_position+length_read)])
# all_data = [change_val(first_df, label_as_int=give_int(genome[initial_position+t].upper())) for t in range(length_read)]
# all_data = [change_val( np.array(create_random_df()), label_as_int=give_int(genome[initial_position+t].upper())) for t in range(length_read)]
all_data = [change_val( np.roll(np.roll(first_df, random.randrange(1, 3), axis=0), random.randrange(1, 3), axis=1),\
 label_as_int=give_int(genome[initial_position+t].upper())) for t in range(length_read)]

print(len(all_data))

total_number_of_frames = len(all_data)



def animate(frame):
    """
    Animation function. Takes the current frame number (to select the potion of
    data to plot) and a line object to update.
    """

    # Not strictly neccessary, just so we know we are stealing these from
    # the global scope
    global all_data, image

    # We want up-to and _including_ the frame'th element
    image.set_array(all_data[frame])

    return image


# Now we can do the plotting!
fig, ax = plt.subplots(1, figsize=(1, 1))
# Remove a bunch of stuff to make sure we only 'see' the actual imshow
# Stretch to fit the whole plane
fig.subplots_adjust(0, 0, 1, 1)
# Remove bounding line
ax.axis("off")

# Initialise our plot. Make sure you set vmin and vmax!
image = ax.imshow(all_data[0], vmin=0, vmax=4, cmap=cmap)
# plt.show()

animation = FuncAnimation(
    # Your Matplotlib Figure object
    fig,
    # The function that does the updating of the Figure
    animate,
    # Frame information (here just frame number)
    np.arange(total_number_of_frames),
    # Extra arguments to the animate function
    fargs=[],
    # Frame-time in ms; i.e. for a given frame-rate x, 1000/x
    interval=1000 / 25
)

# Try to set the DPI to the actual number of pixels you're plotting
animation.save("out_2dgrid.gif")#, dpi=512)