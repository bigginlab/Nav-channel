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
     0.988340676,    0.091031685,    0.121940859,\
     0.096330322,    0.246051773,   -0.964452803,\
    -0.117799774,    0.964967728,    0.234416038,\
    -0.002954572,   -0.000734922, -120.115859985,\
   111.997718811,   93.468566895,   52.190967560,\
   -43.240329742,  283.538330078,   20.000000000 )


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

#sele s1, resi 131-147
#set transparency, 0.5, s1
#show cartoon, s1
#color white, s1

sele vsd, resi 131-232
show cartoon, vsd
color white, vsd

sele pos_charges, vsd and (resn ARG or resn LYS) and not (name C or name N or name O)
show sticks, pos_charges #and ele h and neighbor (ele n+o)
color cb_blue, pos_charges

sele neg_charges, vsd and (resn GLU or resn ASP) and not (name C or name N or name O)
show sticks, neg_charges #and ele h and neighbor (ele n+o)
color cb_red, neg_charges

set cartoon_transparency, 0.5, resi 131-148 # hide s1 helix

# Ensure hetero atoms keep their colours
util.cnc pos_charges
util.cnc neg_charges

# align frames to specific region
intra_fit vsd, 2

# Set resolution and save image
frame 2
ray 2400,2400
png ./VSD_overview_WT.png
