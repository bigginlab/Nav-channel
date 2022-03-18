from sys import argv

my_argv = argv[1:]

# usage: pymol -qt script.pml -- system.pdb 

# load in structure, need to use cmd.load() for arguments
cmd.load(my_argv[0], "system")
cmd.load_traj(my_argv[1], "system") 


# Set some nice CB friendly colours 

set_color cb_orange, [0.8706, 0.5608, 0.0196]
set_color cb_purple, [0.66, 0.35, 0.63]
set_color cb_light_blue, [0.52, 0.75, 0.98]
set_color cb_green, [0.0078, 0.6196, 0.4510]
set_color cb_blue, [0.0039, 0.4510, 0.6980]
set_color cb_red, [0.8353, 0.3686, 0]

set_view (\
     0.873771071,   -0.034619309,    0.485096306,\
     0.449759454,    0.436998397,   -0.778938830,\
    -0.185020015,    0.898795784,    0.397408813,\
    -0.001072153,   -0.000038883,  -72.406684875,\
   110.372161865,   93.032722473,   66.532112122,\
    37.449394226,  107.330337524,  -20.000000000 )

# Start empty and add things
hide everything

# Can either remove all or the non-polar hydrogens
#remove hydrogens
remove (h. and (e. c extend 1))

# Balls and sticks looks better, but some hate it
#set stick_ball, on
#set stick_ball_ratio, 1.5

# Add filter (ambient) and other misc settings
set ambient, 0.4
bg_colour white
#set ray_trace_background, off
set antialias = 1
set ortho = 1
set sphere_mode, 5
set ray_trace_mode, 1
set cartoon_fancy_helices, 1
set cartoon_highlight_color = grey30

sele s3_s4, resi 192-232
show cartoon, s3_s4
color white, s3_s4

sele loop, resi 203-219
hide sticks, loop and name N or name C or name O

sele loop2, resi 209-218 # looking at the residues between S3 and S4
show sticks, loop2
hide sticks, loop2 and name N or name C or name O

sele d211, resi 211
show sticks, d211
hide sticks, d211 and name N or name C or name O
color cb_blue, d211

sele r219_e208, resi 219 or resid 208
show sticks, r219_e208
hide sticks, r219_e208 and name N or name C or name O
color cyan, r219_e208

# Ensure hetero atoms keep their colours
util.cnc loop

# align frames to specific region
intra_fit s3_s4, 50

# Set resolution and save image
frame 1164
ray 2400,2400
png ./D211_no_distortion_WT.png
