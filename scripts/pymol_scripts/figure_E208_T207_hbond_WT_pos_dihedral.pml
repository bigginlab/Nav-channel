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

#set_view (\
#     0.964177728,    0.039527968,    0.262298316,\
#     0.218639717,    0.441484660,   -0.870222986,\
#    -0.150199428,    0.896397412,    0.417027295,\
#     0.000000000,    0.000000000,  -57.191360474,\
#   110.088455200,   88.587799072,   66.027626038,\
#    22.250885010,   92.131813049,  -20.000000000 )

set_view (\
     0.964177728,    0.039527968,    0.262298316,\
     0.218639717,    0.441484660,   -0.870222986,\
    -0.150199428,    0.896397412,    0.417027295,\
     0.000000000,    0.000000000,  -70.888175964,\
   110.088455200,   88.587799072,   66.027626038,\
    35.947689056,  105.828620911,   20.000000000 )


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
show sticks, loop and ele h and neighbor (ele n+o)

sele t207, resi 207
show sticks, t207
hide sticks, t207 and name N or name C or name O
#show sticks, t207 and ele h and neighbor (ele n+o) # show polar hydrogens
color cb_blue, t207

sele r219_e208, resid 219 or resid 208
show sticks, r219_e208
hide sticks, r219_e208 and name N or name C or name O
color cyan, r219_e208

# Ensure hetero atoms keep their colours
util.cnc loop

# align frames to specific region
intra_fit s3_s4, 50

# Set resolution and save image
frame 1079
ray 2400,2400
png ./E208_T207_hbond_WT_pos_dihedral.png
