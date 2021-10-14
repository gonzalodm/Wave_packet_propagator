import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import imageio
import os
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import matplotlib
matplotlib.rcParams['axes.linewidth'] = 2

def hacer_video(cant_fotos):
    dir_name = "images"
    lista_fotos=[] #aca voy a ir guardando las fotos
    for i in range (cant_fotos):
        file_name = os.path.join(dir_name, "out{:05}.png".format(i))
        lista_fotos.append(imageio.imread(file_name))
        print("{} de {} fotos leidas".format(i, cant_fotos-1)) #contador porque tarda un poco

    video_name = os.path.join(dir_name, "animation5.mp4")
    imageio.mimsave(video_name, lista_fotos)
    # imageio.mimsave(video_name, lista_fotos, fps=20)
    print('Video Guardado')
    return

def print_wp(x_grid, istep ,dens_grid, dir_name):

    file_name = os.path.join(dir_name, "out{:05}.png".format(istep))
    print("step",istep)

    fig_dims = (8,6)
    x_aux = [0.4,0.4,0.41,0.41]
    V_aux = [0,25,25,0]
    fig, fig1 = plt.subplots(dpi=120, figsize=fig_dims)
    fig1.plot(x_grid, dens_grid,linewidth=1,linestyle="-")
    # fig1.plot(x_aux, V_aux,linewidth=3.5,linestyle="-",color="k")
    fig1.set_xlim(x_grid[0],x_grid[len(x_grid)-1])
    fig1.set_ylim(-5.0,5.0)
    # fig1.set_ylim(dens_grid[0], dens_grid[len(x_grid)-1])
    # fig1.set_xlim(0.0,30.0)

    fig1.set_xticks([])
    fig1.set_yticks([])
    plt.tight_layout()
    plt.savefig(file_name)


    return

def main():

    files1 = ["time.out", "xgrid.out"]
    files2 = ["psi_re.out"]
    t_len  = 0
    x_len  = 0
    x_grid = []

    entrada = open("./"+files1[0],"r")
    for val in entrada:
        t_len += 1
    entrada.close()

    entrada = open("./"+files1[1],"r")
    for val in entrada:
        x_grid.append(float(val))
    entrada.close()
    x_len = len(x_grid)

    dir_name = "images"
    if not os.path.exists(dir_name):
            os.mkdir(dir_name)

    entrada = open("./"+files2[0],"r")
    for tt in range(t_len):

        dens_grid=[]
        for ii in range(x_len):
            val = entrada.readline()
            dens_grid.append(float(val))

        print_wp(x_grid, tt ,dens_grid, dir_name)

    entrada.close()

    hacer_video(t_len)

    return

main()
