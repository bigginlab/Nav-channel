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
     0.860981584,    0.030409141,    0.507726848,\
     0.477145612,    0.297462702,   -0.826941252,\
    -0.176178798,    0.954241812,    0.241602108,\
    -0.002087608,   -0.001038249, -163.398727417,\
   110.933326721,   90.406120300,   52.317306519,\
  -936.189697266, 1262.999389648,   20.000000000 )

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

sele vsd_not_s3_s4, resi 131-208 or resi 220-232
show cartoon, vsd_not_s3_s4
color white, vsd_not_s3_s4

sele s3_s4, resi 209-219
show cartoon, s3_s4
color violet, s3_s4

#set cartoon_transparency, 0.25, vsd_not_s3_s4 # highlight s3_s4

# align frames to specific region
intra_fit vsd_not_s3_s4, 2

# Set resolution and save image
frame 2
ray 2400,2400
png ./S3_S4_highlight.png
