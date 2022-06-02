# -*- coding:utf-8 -*-
import sys
# sys.path.append('/home/lhaaso/caowy/.mylib/miniconda/miniconda3/envs/myroot/lib/python3.7/site-packages')
sys.path.append('/home/lhaaso/caowy/.mylib/miniconda/miniconda3/')
print('load pakege:')
from math import ceil, sin, cos, radians
from matplotlib import markers
from tqdm import tqdm,trange
pbar = tqdm(total=100)

# import ROOT as rt
from ROOT import TFile
pbar.update(25)
from root_numpy import tree2array
pbar.update(25)
import numpy as np
import pandas as pd

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
# matplotlib.use('TkAgg')
import scipy.ndimage as ndi
pbar.update(25)

from tkinter import *
from tkinter.scrolledtext import *
from tkinter.ttk import *

import io
import os
import glob
import time
import threading
import traceback
pbar.update(25)
pbar.close()

def flatten(items):
    "flatten anything"
    for item in items:
        if isinstance(item, (list, np.ndarray, tuple, np.void)): #any class if you want
            yield from flatten(item)
        else:
            yield item

def repadg():
    global padg
    if ifgd:
        padg = int(np.random.rand()*len(g['x']))
    print('repadg')
    BTt.insert(INSERT,"Repadg:done"+"\n");BTt.see(END)
    return 0

def repadh():
    global padh
    if ifhd:
        padh = int(np.random.rand()*len(h['x']))
    print('repadh')
    BTt.insert(INSERT,"Repadh:done"+"\n");BTt.see(END)
    return 0 

def direction1(phi, theta):
    dx = sin(radians(theta))*cos(radians(phi))*120/5
    dy = sin(radians(theta))*sin(radians(phi))*52/5
    return dx,dy

def direction2(phi, theta):
    dx = sin(radians(theta))*cos(radians(phi))*(detx.max()-detx.min())/5
    dy = sin(radians(theta))*sin(radians(phi))*(dety.max()-dety.min())/5
    return dx,dy

def read_detector():
    global detx
    global dety
    if os.path.exists("./det.npz"):
        det = np.load("./det.npz")
        detx = det['detx']
        dety = det['dety']
    else:
        detx = [e for e in tqdm(flatten(g['x']))]
        dety = [e for e in tqdm(flatten(g['y']))]
        l1['text']='Wait a moment:'
        root.update() 
        detx = np.unique(np.array(detx))
        dety = np.unique(np.array(dety))
        dety=dety[::-1]
        np.savez("det.npz", detx = detx, dety = dety)
    BTt.insert(INSERT,"load detector: Done!!!"+"\n");BTt.see(END)
    return 0

def xy2ij():
    try:
        global gxi
        global gyj
        global hxi
        global hyj
        global padg
        global padh
        
        gxi = []; gyj = []; hxi = []; hyj = []
        for i in tqdm(range(loadcache)):
            ii = i
            if ifgg:
                gxi.append([np.where(detx == x)[0][0] for x in g['x'][ii+padg][0]])
                gyj.append([np.where(dety == y)[0][0] for y in g['y'][ii+padg][0]])
            if ifhh:
                hxi.append([np.where(detx == x)[0][0] for x in h['x'][ii+padh][0]])
                hyj.append([np.where(dety == y)[0][0] for y in h['y'][ii+padh][0]])
        BTt.insert(INSERT,"xy2ij: Done!!!"+"\n");BTt.see(END)
    except Exception as e:
        l1['text']='ERROR!!'
        BTt.insert(INSERT,"ERROR:"+"\n");BTt.see(END)
        BTt.insert(INSERT,str(e.args)+"\n");BTt.see(END)
        BTt.insert(INSERT,"============="+"\n");BTt.see(END)
        BTt.insert(INSERT,traceback.format_exc()+"\n");BTt.see(END)

def rxy2ij(x,y):
    i = int((x-detx.min())/(detx.max()-detx.min())*120)
    j = int((y-dety.min())/(dety.max()-dety.min())*52)
    j = 52-j
    return i,j

def plot():
    global bbj
    global ifg
    global ifh
    global ifgd
    global ifhd
    global cg
    global ch
    global scc
    global ieventg
    global ieventh
    global ieventgg
    global ieventhh
    global padh
    global padg
    fig.clear()
    cm = plt.cm.get_cmap(color.get())
    ifgd = ifg.get()
    ifhd = ifh.get()
    

    if (scc==sc.get()) & (bbj==pltobd.get()):
        ifreread = 0
    else:
        ifreread = 1
        scc = sc.get()
        bbj = pltobd.get()

    sccg='True'
    scch='True'
    scccg = scc
    sccch = scc
    if ifreread:
        try:
            cg[bbj] = np.zeros((loadcache,52,120))+1
            ch[bbj] = np.zeros((loadcache,52,120))+1
            ieventg = []
            ieventh = []
            l1['text']="Searching for events:"
            root.update()
            for item in list(set(Parameters) ^ set(Parametersl)):
                if CV[item].get():
                    if ' '+item+' ' in scc:
                        if ifgd:
                            scccg = scccg.replace(item, ("g['%s'][i+padg][0]" %item))
                            sccg = scccg
                        if ifhd:
                            sccch = sccch.replace(item, ("h['%s'][i+padh][0]" %item))
                            scch = sccch

            while ((not ieventg) and ifgd) or ((not ieventh) and ifhd):
                p1['value']+=1
                root.update()
                if not ieventg:
                    repadg()
                if not ieventh:
                    repadh()

                l1['text']='Searching for events:'
                t2 = threading.Thread(target=xy2ij)
                t2.setDaemon(True)
                t2.start()
                root.update()
                while t2.isAlive():
                    p1['value'] += 20
                    root.update()
                    time.sleep(0.7)

                for i in tqdm(range(loadcache)):
                    if ifgd:
                        if eval(sccg):
                            ieventg.append(i)
                            for j in range(len(gxi[i])):
                                cg[bbj][i][gyj[i][j]][gxi[i][j]] += g[bbj][i+padg][0][j]
                            # cgmin = cg[bbj][i].min()
                            # for ii in range(52):
                            #     for jj in range(120):
                            #         if cg[bbj][i][ii][jj] == 1:
                            #             cg[bbj][i][ii][jj] = cgmin
                    if ifhd:
                        if eval(scch):
                            ieventh.append(i)
                            for j in range(len(hxi[i])):
                                ch[bbj][i][hyj[i][j]][hxi[i][j]] += h[bbj][i+padh][0][j]
                            # chmin = ch[bbj][i].min()
                            # for ii in range(52):
                            #     for jj in range(120):
                            #         if ch[bbj][i][ii][jj] == 1:
                            #             ch[bbj][i][ii][jj] = chmin
        except Exception as e:
            l1['text']='ERROR!!'
            BTt.insert(INSERT,"ERROR:"+"\n");BTt.see(END)
            BTt.insert(INSERT,str(e.args)+"\n");BTt.see(END)
            BTt.insert(INSERT,"============="+"\n");BTt.see(END)
            BTt.insert(INSERT,traceback.format_exc()+"\n");BTt.see(END)

    #Event id
    if ifgd:
        ieventgg = np.random.choice(ieventg)
    if ifhd:
        ieventhh = np.random.choice(ieventh)
    #Smoth sigma, 0 means no smooth
    sigma = smooth.get()

    try:
        if (ifgd and ifhd):
            ax = fig.add_subplot(121)
            bx = fig.add_subplot(122, sharey = ax)
            bx.axes.yaxis.set_visible(False)
            vmin = min(cg[bbj][ieventgg].min()+abs(cg[bbj][ieventgg].min()/5), ch[bbj][ieventhh].min()+abs(ch[bbj][ieventhh].min()/5))
            vmax = max(cg[bbj][ieventgg].max(), ch[bbj][ieventhh].max())
            if bbj == "npe" or bbj == "npe_c":
                norm = colors.LogNorm(clip=True,vmin=vmax/500,vmax=vmax)
            else:
                norm = colors.Normalize(vmin=vmin, vmax=vmax)
            if mode.get() == "image":
                ax.set_xticks(np.linspace(0,120,5))
                ax.set_xticklabels(np.linspace(int(detx.min()),int(detx.max()),5),rotation = 30,fontsize = 'small')
                ax.set_yticks(np.linspace(0,52,5))
                ax.set_yticklabels(np.linspace(int(dety.max()),int(dety.min()),5),rotation = 30,fontsize = 'small')

                bx.set_xticks(np.linspace(0,120,5))
                bx.set_xticklabels(np.linspace(int(detx.min()),int(detx.max()),5),rotation = 30,fontsize = 'small')
                bx.set_yticks(np.linspace(0,52,5))
                bx.set_yticklabels(np.linspace(int(dety.max()),int(dety.min()),5),rotation = 30,fontsize = 'small')

                #g
                a = ax.imshow(ndi.gaussian_filter(cg[bbj][ieventgg], sigma=float(sigma)),aspect='auto',cmap=cm,interpolation='none',norm=norm)

                xc,yc = rxy2ij(g['xc'][ieventgg+padg][0],g['yc'][ieventgg+padg][0])
                xcmc,ycmc = rxy2ij(g['xcmc'][ieventgg+padg][0],g['ycmc'][ieventgg+padg][0])
                dx,dy = direction1(g['phi'][ieventgg+padg][0],g['theta'][ieventgg+padg][0])
                dxmc,dymc = direction1(g['azimc'][ieventgg+padg][0],g['zenmc'][ieventgg+padg][0])
                dx*-1;dy*=-1;dymc*=-1
                ax.scatter(xcmc, ycmc, s=50, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                ax.arrow(xcmc,ycmc,dxmc,dymc,head_width=1, facecolor='blue', label="Real direction", alpha=0.5)
                ax.scatter(xc, yc, s=50, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                ax.arrow(xc,yc,dx,dy,head_width=1, facecolor='red', label="Rec direction", alpha=0.5)

                #h
                b = bx.imshow(ndi.gaussian_filter(ch[bbj][ieventhh], sigma=float(sigma)),aspect='auto',cmap=cm,interpolation='none',norm=norm)

                xc,yc = rxy2ij(h['xc'][ieventhh+padh][0],h['yc'][ieventhh+padh][0])
                xcmc,ycmc = rxy2ij(h['xcmc'][ieventhh+padh][0],h['ycmc'][ieventhh+padh][0])
                dx,dy = direction1(h['phi'][ieventhh+padh][0],h['theta'][ieventhh+padh][0])
                dxmc,dymc = direction1(h['azimc'][ieventhh+padh][0],h['zenmc'][ieventhh+padh][0])
                dx*-1;dy*=-1;dymc*=-1
                bx.scatter(xcmc, ycmc, s=50, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                bx.scatter(xc, yc, s=50, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                bx.arrow(xcmc,ycmc,dxmc,dymc,head_width=1, facecolor='blue', label="Real direction", alpha=0.5)
                bx.arrow(xc,yc,dx,dy,head_width=1, facecolor='red', label="Rec direction", alpha=0.5)
            else:
                #g
                a = ax.scatter(g['x'][ieventgg+padg][0], g['y'][ieventgg+padg][0], s=50, c=g[bbj][ieventgg+padg][0], alpha=0.5, norm=norm, cmap=cm)
                
                dx,dy = direction2(g['phi'][ieventgg+padg][0],g['theta'][ieventgg+padg][0])
                dxmc,dymc = direction2(g['azimc'][ieventgg+padg][0],g['zenmc'][ieventgg+padg][0])

                ax.scatter(g['xcmc'][ieventgg+padg][0], g['ycmc'][ieventgg+padg][0], s=60, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                ax.arrow(g['xcmc'][ieventgg+padg][0],g['ycmc'][ieventgg+padg][0],dxmc,dymc,head_width=5, facecolor='blue', label="Real direction", alpha=0.5)
                ax.scatter(g['xc'][ieventgg+padg][0], g['yc'][ieventgg+padg][0], s=60, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                ax.arrow(g['xc'][ieventgg+padg][0],g['yc'][ieventgg+padg][0],dx,dy,head_width=5, facecolor='red', label="Rec direction", alpha=0.5)

                #h
                b = bx.scatter(h['x'][ieventhh+padh][0], h['y'][ieventhh+padh][0], s=50, c=h[bbj][ieventhh+padh][0], alpha=0.5, norm=norm, cmap=cm)


                dx,dy = direction2(h['phi'][ieventhh+padh][0],h['theta'][ieventhh+padh][0])
                dxmc,dymc = direction2(h['azimc'][ieventhh+padh][0],h['zenmc'][ieventhh+padh][0])

                bx.scatter(h['xcmc'][ieventhh+padh][0], h['ycmc'][ieventhh+padh][0], s=50, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                bx.scatter(h['xc'][ieventhh+padh][0], h['yc'][ieventhh+padh][0], s=50, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
         
                bx.arrow(h['xcmc'][ieventhh+padh][0],h['ycmc'][ieventhh+padh][0],dxmc,dymc,head_width=5, facecolor='blue', label="Real direction", alpha=0.5)
                bx.arrow(h['xc'][ieventhh+padh][0],h['yc'][ieventhh+padh][0],dx,dy,head_width=5, facecolor='red', label="Rec direction", alpha=0.5)

            plt.subplots_adjust(wspace=0, hspace=0)
            ax.set_xlabel("Gamma")
            ax.set_title(("nhit: %d priE: %.2f" %(g['nhit'][ieventgg+padg][0], g['priE'][ieventgg+padg][0])),fontsize=10)
            bx.set_xlabel("Hadorn")
            bx.set_title(("nhit: %d priE: %.2f" %(h['nhit'][ieventhh+padh][0], h['priE'][ieventhh+padh][0])),fontsize=10)
            fig.colorbar(a, ax=[ax, bx], shrink=0.7) #,orientation='horizontal', pad=0.18
        #g
        elif ifgd:
            if mode.get() == "image":
                vmin = cg[bbj][ieventgg].min()+abs(cg[bbj][ieventgg].min()/5)
                vmax = cg[bbj][ieventgg].max()
                if bbj == "npe" or bbj == "npe_c":
                    norm = colors.LogNorm(clip=True,vmin=vmax/500,vmax=vmax)
                else:
                    norm = colors.Normalize(vmin=vmin, vmax=vmax)

                plt.xticks(np.linspace(0,120,5),np.linspace(int(detx.min()),int(detx.max()),5))
                plt.yticks(np.linspace(0,52,5),np.linspace(int(dety.max()),int(dety.min()),5))

                #g
                a = plt.imshow(ndi.gaussian_filter(cg[bbj][ieventgg], sigma=float(sigma)),aspect='auto',cmap=cm,interpolation='none',norm=norm)

                xc,yc = rxy2ij(g['xc'][ieventgg+padg][0],g['yc'][ieventgg+padg][0])
                xcmc,ycmc = rxy2ij(g['xcmc'][ieventgg+padg][0],g['ycmc'][ieventgg+padg][0])
                dx,dy = direction1(g['phi'][ieventgg+padg][0],g['theta'][ieventgg+padg][0])
                dxmc,dymc = direction1(g['azimc'][ieventgg+padg][0],g['zenmc'][ieventgg+padg][0])
                dx*-1;dy*=-1;dymc*=-1
                plt.scatter(xcmc, ycmc, s=50, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                plt.arrow(xcmc,ycmc,dxmc,dymc,head_width=1, facecolor='blue', label="Real direction", alpha=0.5)
                plt.scatter(xc, yc, s=50, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                plt.arrow(xc,yc,dx,dy,head_width=1, facecolor='red', label="Rec direction", alpha=0.5)
            else:
                vmin = g[bbj][ieventgg+padg][0].min()+abs(cg[bbj][ieventgg].min()/5)
                vmax = g[bbj][ieventgg+padg][0].max()
                if bbj == "npe" or bbj == "npe_c":
                    norm = colors.LogNorm(clip=True,vmin=vmax/500,vmax=vmax)
                else:
                    norm = colors.Normalize(vmin=vmin, vmax=vmax)
                a = plt.scatter(g['x'][ieventgg+padg][0], g['y'][ieventgg+padg][0], s=50, c=g[bbj][ieventgg+padg][0], alpha=0.5, norm=norm, cmap=cm)

                dx,dy = direction2(g['phi'][ieventgg+padg][0],g['theta'][ieventgg+padg][0])
                dxmc,dymc = direction2(g['azimc'][ieventgg+padg][0],g['zenmc'][ieventgg+padg][0])

                plt.scatter(g['xcmc'][ieventgg+padg][0], g['ycmc'][ieventgg+padg][0], s=60, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                plt.arrow(g['xcmc'][ieventgg+padg][0],g['ycmc'][ieventgg+padg][0],dxmc,dymc,head_width=5, facecolor='blue', label="Real direction", alpha=0.5)
                plt.scatter(g['xc'][ieventgg+padg][0], g['yc'][ieventgg+padg][0], s=60, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                plt.arrow(g['xc'][ieventgg+padg][0],g['yc'][ieventgg+padg][0],dx,dy,head_width=5, facecolor='red', label="Rec direction", alpha=0.5)
            plt.xlabel("Gamma")
            plt.title(("nhit: %d priE: %.2f" %(g['nhit'][ieventgg+padg][0], g['priE'][ieventgg+padg][0])),fontsize=10)
            fig.colorbar(a, shrink=0.7) #, orientation='horizontal', pad=0.18
        #h
        elif ifhd:
            if mode.get() == "image":
                vmin = ch[bbj][ieventhh].min()+abs(ch[bbj][ieventhh].min()/5)
                vmax = ch[bbj][ieventhh].max()
                if bbj == "npe" or bbj == "npe_c":
                    norm = colors.LogNorm(clip=True,vmin=vmax/500,vmax=vmax)
                else:
                    norm = colors.Normalize(vmin=vmin, vmax=vmax)
                plt.xticks(np.linspace(0,120,5),np.linspace(int(detx.min()),int(detx.max()),5))
                plt.yticks(np.linspace(0,52,5),np.linspace(int(dety.max()),int(dety.min()),5))
                b = plt.imshow(ndi.gaussian_filter(ch[bbj][ieventhh], sigma=float(sigma)),aspect='auto',cmap=cm,interpolation='none',norm=norm)

                xc,yc = rxy2ij(h['xc'][ieventhh+padh][0],h['yc'][ieventhh+padh][0])
                xcmc,ycmc = rxy2ij(h['xcmc'][ieventhh+padh][0],h['ycmc'][ieventhh+padh][0])
                dx,dy = direction1(h['phi'][ieventhh+padh][0],h['theta'][ieventhh+padh][0])
                dxmc,dymc = direction1(h['azimc'][ieventhh+padh][0],h['zenmc'][ieventhh+padh][0])
                dx*-1;dy*=-1;dymc*=-1
                plt.scatter(xcmc, ycmc, s=50, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                plt.scatter(xc, yc, s=50, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                plt.arrow(xcmc,ycmc,dxmc,dymc,head_width=1, facecolor='blue', label="Real direction", alpha=0.5)
                plt.arrow(xc,yc,dx,dy,head_width=1, facecolor='red', label="Rec direction", alpha=0.5)
            else:
                vmin = h[bbj][ieventhh+padh][0].min()+abs(ch[bbj][ieventhh].min()/5)
                vmax = h[bbj][ieventhh+padh][0].max()
                if bbj == "npe" or bbj == "npe_c":
                    norm = colors.LogNorm(clip=True,vmin=vmax/500,vmax=vmax)
                else:
                    norm = colors.Normalize(vmin=vmin, vmax=vmax)
                b = plt.scatter(h['x'][ieventhh+padh][0], h['y'][ieventhh+padh][0], s=50, c=h[bbj][ieventhh+padh][0], alpha=0.5, norm=norm, cmap=cm)
                dx,dy = direction2(h['phi'][ieventhh+padh][0],h['theta'][ieventhh+padh][0])
                dxmc,dymc = direction2(h['azimc'][ieventhh+padh][0],h['zenmc'][ieventhh+padh][0])

                plt.scatter(h['xcmc'][ieventhh+padh][0], h['ycmc'][ieventhh+padh][0], s=50, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                plt.scatter(h['xc'][ieventhh+padh][0], h['yc'][ieventhh+padh][0], s=50, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
         
                plt.arrow(h['xcmc'][ieventhh+padh][0],h['ycmc'][ieventhh+padh][0],dxmc,dymc,head_width=5, facecolor='blue', label="Real direction", alpha=0.5)
                plt.arrow(h['xc'][ieventhh+padh][0],h['yc'][ieventhh+padh][0],dx,dy,head_width=5, facecolor='red', label="Rec direction", alpha=0.5)
            plt.xlabel("Hadorn")
            plt.title(("nhit: %d priE: %.2f" %(h['nhit'][ieventhh+padh][0], h['priE'][ieventhh+padh][0])),fontsize=10)
            fig.colorbar(b, shrink=0.7)
        else:
            BTt.insert(INSERT,"Error: \n");BTt.see(END)
            BTt.insert(INSERT,"Please choose a file.\n");BTt.see(END)
        
        titg = ""
        tith = ""
        for item in Parameterso:
            if CV[item].get():
                if item == "pinc":
                    if ifgd:
                        titg += str(item)+(("2:%.2f "%g[str(item)][ieventgg+padg][0][2]) + "\t")
                    if ifhd:
                        tith += str(item)+(("2:%.2f "%h[str(item)][ieventhh+padh][0][2]) +"\t")
                elif item == "nhit":
                    if ifgd:
                        titg += str(item)+((":%d "%g[str(item)][ieventgg+padg][0]) + "\t")
                    if ifhd:
                        tith += str(item)+((":%d " %h[str(item)][ieventhh+padh][0]) + "\t")
                else:
                    if ifgd:
                        titg += str(item)+((":%.2f " %g[str(item)][ieventgg+padg][0]) +"\t")
                    if ifhd:
                        tith += str(item)+((":%.2f " %h[str(item)][ieventhh+padh][0] )+"\t")
        if ifgd:
            BTt.insert(INSERT,"Event information of Gamma: \n");BTt.see(END)
            BTt.insert(INSERT,titg+"\n");BTt.see(END)
        if ifhd:
            BTt.insert(INSERT,"Event information of Hadorn: \n");BTt.see(END)
            BTt.insert(INSERT,tith+"\n");BTt.see(END)
       
        # plt.contour(cgpe[ievent], [cgpe[ievent].max()-400, cgpe[ievent].max()-300, cgpe[ievent].max()-200])
        plt.legend(loc="upper left", bbox_to_anchor=(0.9,1.16),framealpha=0.3)
        canvas_spice.draw()
        toolbar.update()
        plt.ion()
        l1['text']='Drawing completed in plot! '
        BTt.insert(INSERT,"Drawing completed in plot!"+"\n");BTt.see(END)
        root.update()
    except Exception as e:
        l1['text']='ERROR!!'
        BTt.insert(INSERT,"ERROR:"+"\n");BTt.see(END)
        BTt.insert(INSERT,str(e.args)+"\n");BTt.see(END)
        BTt.insert(INSERT,"============="+"\n");BTt.see(END)
        BTt.insert(INSERT,traceback.format_exc()+"\n");BTt.see(END)
    return 0

def load():
    global l1
    global p1
    global pltobd
    global Parameters
    global ifgg
    global ifhh
    global fgname
    global fhname
    global tgamma
    global thadron
    global dg
    global dh

    l1['text']="Load data:"
    BTt.insert(INSERT,"Load data:"+"\n");BTt.see(END)
    
    pltobj.destroy()
    pltobd = Combobox(frame4,width=15,values=[i for i in Parametersll if CV[i].get()])
    pltobd.grid(row=3, column=0)
    root.update()

    fgname = e1.get()
    fhname = e2.get()

    ifgg = not(fgname == "")
    ifhh = not(fhname == "")

    ifg.set(ifgg)
    ifh.set(ifhh)

    ggg={}
    hhh={}

    try:
        if not (fgname.endswith('.npz') or fhname.endswith('.npz')):
            if ifgg:
                fgamma=TFile(fgname)
                tgamma = fgamma.Get("rec")
                Parameters= list(set(Parameters).union(set([str(i).split()[1] for i in list(tgamma.GetListOfBranches())])))

            if ifhh:
                fhadron=TFile(fhname)
                thadron = fhadron.Get("rec")
                Parameters= list(set(Parameters).union(set([str(i).split()[1] for i in list(thadron.GetListOfBranches())])))
            for item in C.keys():
                if CV[item].get():
                    print("Load "+item)
                    if ifgg:
                        l1['text']="Load "+item+" of gamma:"
                        root.update()
                        ggg[item]=threading.Thread(target=gload,args=(item,))
                        ggg[item].setDaemon(True)
                        ggg[item].start()
                    if ifhh:
                        l1['text']="Load "+item+" of hadron:"
                        root.update()
                        hhh[item]=threading.Thread(target=hload,args=(item,))
                        hhh[item].setDaemon(True)
                        hhh[item].start()
                        
            p1['value']=100
            l1['text']='Done!!!'
            root.update()
            BTt.insert(INSERT,"Load data: Done!!!"+"\n");BTt.see(END)
            time.sleep(1)
        else:
            for item in C.keys():
                if CV[item].get():
                    if ifgg:
                        dg = np.load(fgname,allow_pickle=True)
                        print("Load "+item)
                        l1['text']="Load "+item+" of gamma:"
                        root.update()
                        ggg[item]=threading.Thread(target=dgload, args=(item,))
                        ggg[item].setDaemon(True)
                        ggg[item].start()

                    if ifhh:
                        dh = np.load(fhname,allow_pickle=True) 
                        l1['text']="Load "+item+" of hadron:"
                        root.update()
                        hhh[item]=threading.Thread(target=dhload, args=(item,))
                        hhh[item].setDaemon(True)
                        hhh[item].start()

                p1['value']=100
                l1['text']='Done!!!'
                root.update()
                BTt.insert(INSERT,"Load data: Done!!!"+"\n");BTt.see(END)
    except Exception as e:
        l1['text']='ERROR!!'
        BTt.insert(INSERT,"ERROR:"+"\n");BTt.see(END)
        BTt.insert(INSERT,str(e.args)+"\n");BTt.see(END)
        BTt.insert(INSERT,"============="+"\n");BTt.see(END)
        BTt.insert(INSERT,traceback.format_exc()+"\n");BTt.see(END)
        
    l1['text']='Load detecter:'
    root.update() 

    t = threading.Thread(target=read_detector)
    t.setDaemon(True)
    t.start()
    while t.isAlive():
        p1['value'] += 20
        root.update()
        time.sleep(0.7)

    p1['value']=100
    l1['text']='Done!!!'
    root.update()

    p1['value']=100
    l1['text']='Now you can plot.'
    root.update()
    return 0

def gload(item):
    g['%s'%item] = tree2array(tgamma, branches=['%s' %item])
    p1['value']+=30
    root.update()
    BTt.insert(INSERT,"Load "+item+" of gamma: Done!!!"+"\n");BTt.see(END)

def hload(item):
    h['%s'%item] = tree2array(thadron, branches=['%s' %item])
    p1['value']+=30
    root.update()
    BTt.insert(INSERT,"Load "+item+" of hadron: Done!!!"+"\n");BTt.see(END)

def dgload(item):
    g['%s'%item] = dg['%s'%item]
    p1['value']+=30
    BTt.insert(INSERT,"Load "+item+" of gamma: Done!!!"+"\n");BTt.see(END)
    root.update()

def dhload(item):
    h['%s'%item] = dh['%s'%item]
    p1['value']+=30
    BTt.insert(INSERT,"Load "+item+" of hadron: Done!!!"+"\n");BTt.see(END)
    root.update()
    
if __name__ == "__main__":
    # Parameters = ['x', 'y', 'npe', 't', 'tr', 'dr', 'cpx', 'pinc', 'priE', 'nhit']
    Parameters = ['x','y','npe','t','xcmc','ycmc','irun','evt','nhit','nfit','sump1','sump2','mjd','sigma','theta','phi','xc','yc','ra','dec','fee','ch','iconf','npe_c','dr','tr','cpx','pinc','csr','dcore','dangle','priE','zenmc','azimc','h1','weight','dangle_real','dcore_real']
    Parametersl = ['x', 'y', 'fee', 'ch', 'iconf', 'npe', 'npe_c', 'dr', 'tr', 't']
    Parametersll = ['npe', 'npe_c', 'tr', 't', 'dr', 'fee', 'ch', 'iconf']
    Parametersm = ['x','y','xcmc','ycmc','xc','yc','zenmc','azimc','theta','phi',"nhit",'priE']
    Parametersu = ['npe', 'npe_c', 't', 'tr', 'dr', 'fee', 'ch', 'iconf', 'irun', 'evt', 'nfit', 'sump1', 'sump2', 'cpx', 'pinc',  'mjd', 'dcore', 'dangle','dcore_real', 'dangle_real', 'ra', 'dec', 'sigma', 'csr', 'h1', 'weight']
    Parameterso = ['theta','phi','zenmc','azimc','irun', 'evt', 'nfit', 'sump1', 'sump2', 'cpx', 'pinc',  'mjd', 'dcore', 'dangle','dcore_real', 'dangle_real', 'ra', 'dec', 'sigma', 'csr', 'h1', 'weight']
    # pl = len(list(set(Parameters) ^ set(Parametersm)))
    pl = len(Parametersu)
    global scc
    scc = ""

    loadcache = 1000
    bbj=''
    g = {}
    h = {}
    cg = {}
    ch = {}
    ieventg = []
    ieventh = []
    ieventgg=0
    ieventhh=0


    root = Tk()
    ww=500
    hh=450
    root.title('WCDA-Viewer V0.3.15')
    # 获取屏幕 宽、高
    ws = root.winfo_screenwidth()
    hs = root.winfo_screenheight()
    # 计算 x, y 位置
    xx = (ws/2) - (ww/2)
    yy = (hs/2) - (hh/2)
    root.geometry('%dx%d+%d+%d' % (ww, hh, xx, yy))
    # root.geometry("500x450")
    root.resizable(0,0)


    s=Style()
    s.theme_names()
    s.theme_use('clam')

    notebook = Notebook(root,height=350, width=500)  # notebook属于ttk

    frame1 = Frame(notebook,height=300,width=500)

    Label(frame1, text = '1. Load Data', font = ('microsoft yahei', 16, 'bold')).grid(row=0, column=0, columnspan=2)

    Label(frame1, text='The ROOT/npz file of Gamma:   ').grid(row=1, column=0, sticky=W)

    e1 = Combobox(frame1, values=['/eos/user/z/zwang/wcda/data/mc/g4s2_rec/gamma_rec_143800t.root', 'datag.npz'], width=37)
    e1.grid(row=1, column=1)
    Label(frame1, text='The ROOT/npz file of Hadron:  ').grid(row=2, column=0, sticky=W)
    e2 = Combobox(frame1, values=['/eos/user/z/zwang/wcda/data/mc/g4s2_rec/hadron_rec_1970s.root', 'datah.npz'], width=37)
    e2.grid(row=2, column=1)

    C={}
    CV={}
    for p in Parametersm:
        CV[p]=IntVar(value=1)
        C[p]=0
    canvas=Canvas(frame1,width=500,height=60,scrollregion=(0,0,1600,200))
    canvas.grid(row=3, column=0, columnspan=2)
    frame11 = Frame(canvas, width=500, height=60)
    frame11.place()
    hbar=Scrollbar(canvas,orient=HORIZONTAL)
    hbar.place(x =0,y=50,width=500,height=10)
    hbar.configure(command=canvas.xview)
    canvas.config(xscrollcommand=hbar.set)
    i=0
    for it in Parametersu:
        CV['%s'%it] = IntVar(value=0)
        C['%s'%it] = Checkbutton(frame11, text = ("%s" %it), onvalue = 1, offvalue = 0, width=13, variable=CV['%s'%it])
        C['%s'%it].grid(row=(0 if i < ceil(pl/2) else 1),column=(i if i < ceil(pl/2) else i-ceil(pl/2)), sticky = W)
        i+=1
    canvas.create_window((0,2), window=frame11,anchor=NW)

    button = Button(frame1, text='Load Data', command=load)
    button.grid(row=4, columnspan=2)

    Label(frame1, text = '2. Plot Data', font = ('microsoft yahei', 16, 'bold')).grid(row=5, column=0, columnspan=2)
    frame4 = Frame(frame1,height=300,width=500)
    frame4.grid(row=6,columnspan=2)
    Label(frame4, text='Mode:').grid(row=0, column=0, sticky=W)
    Label(frame4, text='Color:').grid(row=0, column=1, sticky=W)
    Label(frame4, text='Smooth Radius:').grid(row=0, column=2, sticky=W)
    color = Combobox(frame4,values=list(plt.cm.datad.keys())) #["magma","Greens","Blues", "OrRd", "PiYG"]
    color.grid(row=1, column=1)
    mode = Combobox(frame4,values=["image","scatter"])
    mode.grid(row=1, column=0)
    smooth = Spinbox(frame4, from_=0,to=5,increment=0.1)
    smooth.grid(row=1, column=2)
    Label(frame4, text='What obj do you want to plot:').grid(row=2, columnspan=3, sticky=W)
    pltobj = Combobox(frame4,width=15,values=[list(C.keys())[i] for i in range(2,6) if CV[list(C.keys())[i]].get()])
    pltobj.grid(row=3, column=0)
    ifh = IntVar(value=0)
    ifg = IntVar(value=0)
    bg = Checkbutton(frame4, text = "Gamma?", onvalue = 1, offvalue = 0, width=13, variable=ifg)
    bh = Checkbutton(frame4, text = "Hadron?", onvalue = 1, offvalue = 0, width=13, variable=ifh)
    bg.grid(row=3, column=1)
    bh.grid(row=3, column=2)


    global sc
    sc = StringVar(value="( nhit >= 1000) & ( priE > 10)")
    Label(frame4, text='Selection criteria:').grid(row=4, columnspan=3, sticky=W)
    select = Entry(frame4,width=52,textvariable=sc)
    select.grid(row=5, columnspan=3)

    button = Button(frame1, text='Plot Data', command=plot)
    button.grid(row=7, columnspan=2)



    frame2 = Frame(notebook)
    notebook.add(frame1, text="Setting")
    notebook.add(frame2, text="Plot")
    notebook.grid(row=0, column=0, columnspan=2)

    global p1
    global l1
    global BTt
    l1=Label(root, text='Wait for the operation:')
    l1.grid(row=1,column=0, sticky=W)
    p1 = Progressbar(root,length=200,maximum=110,mode="indeterminate",
        orient=HORIZONTAL)
    p1.grid(row=1, column=1, sticky=W)

    BTt = ScrolledText(root,width=67,height=3,relief="solid",bg="#DCDAD5")
    BTt.grid(row=2, columnspan=2, sticky=SW, ipady=5)
    BTt.insert(INSERT,"This is the output log:"+"\n");BTt.see(END)

    #plt fig
    global fig
    global canvas_spice
    global toolbar
    global frame22
    fig = plt.figure(figsize=(6.6,3.95),dpi=80) 

    canvas_spice = FigureCanvasTkAgg(fig,frame2)
    frame22 = Frame(frame2)
    frame22.grid(row=1)
    toolbar = NavigationToolbar2Tk(canvas_spice, frame22)
    canvas_spice.get_tk_widget().grid(row=0,column=0,columnspan=2)
    canvas_spice._tkcanvas.grid(row=0,column=0,columnspan=2)      #.pack(side=TOP, fill=BOTH, expand=1)
    Button(frame2, text='Next>', command=plot).grid(row=1, column=1)


    print('Show windows')

    root.mainloop()
