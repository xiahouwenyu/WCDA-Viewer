# -*- coding:utf-8 -*-
print('load pakege:')
from math import ceil, sin, cos, radians, floor
from matplotlib import markers
from tqdm import tqdm,trange
pbar = tqdm(total=100)

from ROOT import TFile
pbar.update(25)
from root_numpy import tree2array
pbar.update(25)
import numpy as np

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
matplotlib.use('TkAgg')
import scipy.ndimage as ndi
pbar.update(25)

from tkinter import *
from tkinter.scrolledtext import *
from tkinter.ttk import *

import io
import os
import sys
import glob
import time
import datetime
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
    dx = sin(radians(theta))*cos(radians(phi))*60/5
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
            if ifgd:
                gxi.append([floor(np.where(detx == x)[0][0]/2) for x in g['x'][ii+padg][0]])
                gyj.append([np.where(dety == y)[0][0] for y in g['y'][ii+padg][0]])
            if ifhd:
                hxi.append([floor(np.where(detx == x)[0][0]/2) for x in h['x'][ii+padh][0]])
                hyj.append([np.where(dety == y)[0][0] for y in h['y'][ii+padh][0]])
        BTt.insert(INSERT,"xy2ij: Done!!!"+"\n");BTt.see(END)
    except Exception as e:
        l1['text']='ERROR!!'
        BTt.insert(INSERT,"ERROR:"+"\n");BTt.see(END)
        BTt.insert(INSERT,str(e.args)+"\n");BTt.see(END)
        BTt.insert(INSERT,"============="+"\n");BTt.see(END)
        BTt.insert(INSERT,traceback.format_exc()+"\n");BTt.see(END)

def rxy2ij(x,y):
    i = int((x-detx.min())/(detx.max()-detx.min())*60)
    j = int((y-dety.min())/(dety.max()-dety.min())*52)
    j = 52-j
    return i,j

def plot():
    global bbj
    global ifg
    global ifh
    global ifgd
    global ifhd
    global ifrg
    global ifrh
    global cg
    global ch
    global scc
    global ieventg
    global ieventh
    global ieventgg
    global ieventhh
    global padh
    global padg
    global CV
    global newWindow
    global fig2
    plt.close()
    plt.figure(figsize=(6.6,3.95),dpi=80)
    fm = plt.get_current_fig_manager()
    fm.canvas.figure = fig
    fig.canvas = fm.canvas
    fig.clear()
    cm = plt.cm.get_cmap(color.get())
    ifgd = ifg.get()
    ifhd = ifh.get()
    if (not scc==sc.get()) or (not bbj==pltobd.get()):
        ifrg=0
        ifrh=0

    if ((scc==sc.get()) & (bbj==pltobd.get()) and (not (ifgd==1 and ifrg==0)) and (not (ifhd==1 and ifrh==0))):
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
            if ifgd:
                cg[bbj] = np.zeros((loadcache,52,60))+1
                ieventg = []
            if ifhd:
                ch[bbj] = np.zeros((loadcache,52,60))+1
                ieventh = []
            l1['text']="Searching for events:"
            root.update()
            for item in Parameterso:
                if CV[item].get():
                    if ' '+item+' ' in scc:
                        if ifgd:
                            scccg = scccg.replace(item, ("g['%s'][i+padg][0]" %item))
                            # scccg = re.sub(item, ("g['%s'][i+padg][0]" %item), scccg)
                        if ifhd:
                            sccch = sccch.replace(item, ("h['%s'][i+padh][0]" %item))
                            # sccch = re.sub(item, ("h['%s'][i+padh][0]" %item), sccch)
            sccg = scccg
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
                    if ifgd and (not ifrg):
                        if eval(sccg):
                            ieventg.append(i)
                            for j in range(len(gxi[i])):
                                cg[bbj][i][gyj[i][j]][gxi[i][j]] += g[bbj][i+padg][0][j]
                            # cgmin = cg[bbj][i].min()
                            # for ii in range(52):
                            #     for jj in range(120):
                            #         if cg[bbj][i][ii][jj] == 1:
                            #             cg[bbj][i][ii][jj] = cgmin
                    if ifhd and (not ifrh):
                        if eval(scch):
                            ieventh.append(i)
                            for j in range(len(hxi[i])):
                                ch[bbj][i][hyj[i][j]][hxi[i][j]] += h[bbj][i+padh][0][j]
                if ieventg:
                    ifrg=1
                if ieventh:
                    ifrh=1
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

    t0=datetime.datetime(1858,11,17,0,0,0,0)
    try:
        titg = ""
        tith = ""
        for item in Parametersoo:
            if not mc:
                if "mc" in item or item=='priE':
                    continue
            if CV[item].get():
                if item == "pinc":
                    if ifgd:
                        titg += str(item)+(("2:%.2f "%g[str(item)][ieventgg+padg][0][2]) + "\n")
                    if ifhd:
                        tith += str(item)+(("2:%.2f "%h[str(item)][ieventhh+padh][0][2]) +"\n")
                elif item == "mjd":
                    if ifgd:
                        tg = t0+datetime.timedelta(days=g[str(item)][ieventgg+padg][0]+1)
                        tg = datetime.datetime.date(tg)
                        titg += str(item)+((":%s " %tg) +"\n")
                    if ifhd:
                        th = t0+datetime.timedelta(days=h[str(item)][ieventgg+padg][0]+1)
                        th = datetime.datetime.date(th)
                        tith += str(item)+((":%s " %th) +"\n")
                elif item == "evt":
                    if ifgd:
                        titg += str(item)+((":%d " %g[str(item)][ieventgg+padg][0]) +"\n")
                    if ifhd:
                        tith += str(item)+((":%d " %h[str(item)][ieventhh+padh][0] )+"\n")
                else:
                    if ifgd:
                        titg += str(item)+((":%.2f " %g[str(item)][ieventgg+padg][0]) +"\n")
                    if ifhd:
                        tith += str(item)+((":%.2f " %h[str(item)][ieventhh+padh][0] )+"\n")

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
                ax.set_yticks(np.linspace(0,60,5))
                ax.set_yticklabels(np.linspace(int(detx.max()),int(detx.min()),5),rotation = 30,fontsize = 'small')
                ax.set_xticks(np.linspace(0,52,5))
                ax.set_xticklabels(np.linspace(int(dety.max()),int(dety.min()),5),rotation = 30,fontsize = 'small')

                bx.set_yticks(np.linspace(0,60,5))
                bx.set_yticklabels(np.linspace(int(detx.max()),int(detx.min()),5),rotation = 30,fontsize = 'small')
                bx.set_xticks(np.linspace(0,52,5))
                bx.set_xticklabels(np.linspace(int(dety.max()),int(dety.min()),5),rotation = 30,fontsize = 'small')

                rect11 = plt.Rectangle((0,0),22,60,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect21 = plt.Rectangle((23,0),29,29.5,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect31 = plt.Rectangle((23,31),29,29,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3)
                rect12 = plt.Rectangle((0,0),22,60,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect22 = plt.Rectangle((23,0),29,29.5,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect32 = plt.Rectangle((23,31),29,29,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3)

                #g
                a = ax.imshow(ndi.gaussian_filter(cg[bbj][ieventgg].T, sigma=float(sigma)),cmap=cm,interpolation='none',norm=norm,aspect='equal') #
                # limits = ax.axis()

                xc,yc = rxy2ij(g['xc'][ieventgg+padg][0],g['yc'][ieventgg+padg][0])
                dx,dy = direction1(g['phi'][ieventgg+padg][0],g['theta'][ieventgg+padg][0])
                dx*-1;dy*=-1
                if mc:
                    xcmc,ycmc = rxy2ij(g['xcmc'][ieventgg+padg][0],g['ycmc'][ieventgg+padg][0])
                    dxmc,dymc = direction1(g['azimc'][ieventgg+padg][0],g['zenmc'][ieventgg+padg][0])
                    dymc*=-1
                    ax.scatter(ycmc, xcmc, s=50, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                    ax.arrow(ycmc,xcmc,dymc,dxmc,head_width=1, facecolor='blue', label="Real direction", alpha=0.5)
                ax.scatter(yc, xc, s=50, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                ax.arrow(yc,xc,dy,dx,head_width=1, facecolor='red', label="Rec direction", alpha=0.5)
                ax.add_patch(rect11);ax.add_patch(rect21);ax.add_patch(rect31)
                ax.text(0,0,titg,ha='left',va="bottom",wrap=True,family='serif', style='italic',bbox={'facecolor': 'white','edgecolor':'none','alpha': 0.3,'pad': 0.1,'boxstyle':'round4'}, transform=ax.transAxes)
                # if mc:
                #     ax.set_ylim([min([xcmc,xcmc+dxmc,xc,xc+dx,-5]), max([xcmc,xcmc+dxmc,xc,xc+dx,55])])
                #     ax.set_xlim([min([ycmc,ycmc+dymc,yc,yc+dy,-5]), max([ycmc,ycmc+dymc,yc,yc+dy,65])])
                # ax.set_ylim([min([xc,xc+dx,-5]), max([xc,xc+dx,55])])
                # ax.set_xlim([min([yc,yc+dy,-5]), max([yc,yc+dy,65])])

                #h
                b = bx.imshow(ndi.gaussian_filter(ch[bbj][ieventhh].T, sigma=float(sigma)),cmap=cm,interpolation='none',norm=norm,aspect='equal') #
                # limits = bx.axis()

                xc,yc = rxy2ij(h['xc'][ieventhh+padh][0],h['yc'][ieventhh+padh][0])
                dx,dy = direction1(h['phi'][ieventhh+padh][0],h['theta'][ieventhh+padh][0])
                dx*-1;dy*=-1

                if mc:
                    xcmc,ycmc = rxy2ij(h['xcmc'][ieventhh+padh][0],h['ycmc'][ieventhh+padh][0])
                    dxmc,dymc = direction1(h['azimc'][ieventhh+padh][0],h['zenmc'][ieventhh+padh][0])
                    dymc*=-1
                    bx.scatter(ycmc, xcmc, s=50, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                    bx.arrow(ycmc,xcmc,dymc,dxmc,head_width=1, facecolor='blue', label="Real direction", alpha=0.5)

                bx.scatter(yc, xc, s=50, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                bx.arrow(yc,xc,dy,dx,head_width=1, facecolor='red', label="Rec direction", alpha=0.5)
                bx.add_patch(rect12);bx.add_patch(rect22);bx.add_patch(rect32)
                bx.text(1,0,tith,ha='left',va="bottom", wrap=True,family='serif', style='italic',bbox={'facecolor': 'white','edgecolor':'none','alpha': 0.3,'pad': 0.1,'boxstyle':'round4'}, transform=ax.transAxes)
                # bx.axis(limits)
                # if mc:
                #     bx.set_ylim([min([xcmc,xcmc+dxmc,xc,xc+dx,-5]), max([xcmc,xcmc+dxmc,xc,xc+dx,55])])
                #     bx.set_xlim([min([ycmc,ycmc+dymc,yc,yc+dy,-5]), max([ycmc,ycmc+dymc,yc,yc+dy,65])])
                # bx.set_ylim([min([xc,xc+dx,-5]), max([xc,xc+dx,55])])
                # bx.set_xlim([min([yc,yc+dy,-5]), max([yc,yc+dy,65])])
            else:
                rect11 = plt.Rectangle((-135,-150),150,150,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect21 = plt.Rectangle((-135,5),150,150,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect31 = plt.Rectangle((20,-147.5),110,300,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3)

                rect12 = plt.Rectangle((-135,-150),150,150,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect22 = plt.Rectangle((-135,5),150,150,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect32 = plt.Rectangle((20,-147.5),110,300,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3)
                #g
                a = ax.scatter(g['y'][ieventgg+padg][0], g['x'][ieventgg+padg][0], s=50, c=g[bbj][ieventgg+padg][0], alpha=0.5, norm=norm, cmap=cm)
                # limits = ax.axis()
                
                dx,dy = direction2(g['phi'][ieventgg+padg][0],g['theta'][ieventgg+padg][0])

                if mc:
                    dxmc,dymc = direction2(g['azimc'][ieventgg+padg][0],g['zenmc'][ieventgg+padg][0])
                    ax.scatter(g['ycmc'][ieventgg+padg][0], g['xcmc'][ieventgg+padg][0], s=60, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                    ax.arrow(g['ycmc'][ieventgg+padg][0],g['xcmc'][ieventgg+padg][0],dymc,dxmc,head_width=5, facecolor='blue', label="Real direction", alpha=0.5)
                ax.scatter(g['yc'][ieventgg+padg][0], g['xc'][ieventgg+padg][0], s=60, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                ax.arrow(g['yc'][ieventgg+padg][0],g['xc'][ieventgg+padg][0],dy,dx,head_width=5, facecolor='red', label="Rec direction", alpha=0.5)
                ax.add_patch(rect11);ax.add_patch(rect21);ax.add_patch(rect31)
                ax.text(0, 0,titg,ha='left',va="bottom", wrap=True,family='serif', style='italic',bbox={'facecolor': 'white','edgecolor':'none','alpha': 0.3,'pad': 0.1,'boxstyle':'round4'}, transform=ax.transAxes)
                # if mc:
                #     ax.set_ylim([min([g['xcmc'][ieventgg+padg][0],g['xcmc'][ieventgg+padg][0]+dxmc,g['xc'][ieventgg+padg][0],-130]), max([g['xcmc'][ieventgg+padg][0],g['xcmc'][ieventgg+padg][0]+dxmc,g['xc'][ieventgg+padg][0],g['xc'][ieventgg+padg][0]+dx,130])])
                #     ax.set_xlim([min(g['ycmc'][ieventgg+padg][0],g['ycmc'][ieventgg+padg][0]+dymc,g['yc'][ieventgg+padg][0],g['yc'][ieventgg+padg][0]+dy,-150), max(g['ycmc'][ieventgg+padg][0],g['ycmc'][ieventgg+padg][0]+dymc,g['yc'][ieventgg+padg][0],g['yc'][ieventgg+padg][0]+dy,150)])
                # ax.set_ylim([min([g['xc'][ieventgg+padg][0],g['xc'][ieventgg+padg][0]+dx,-130]), max([g['xc'][ieventgg+padg][0],g['xc'][ieventgg+padg][0]+dx,55])])
                # ax.set_xlim([min(g['yc'][ieventgg+padg][0],g['yc'][ieventgg+padg][0]+dy,-150), max(g['yc'][ieventgg+padg][0],g['yc'][ieventgg+padg][0]+dy,150)])

                #h
                b = bx.scatter(h['y'][ieventhh+padh][0], h['x'][ieventhh+padh][0], s=50, c=h[bbj][ieventhh+padh][0], alpha=0.5, norm=norm, cmap=cm)
                # limits = bx.axis()


                dx,dy = direction2(h['phi'][ieventhh+padh][0],h['theta'][ieventhh+padh][0])

                if mc:
                    dxmc,dymc = direction2(h['azimc'][ieventhh+padh][0],h['zenmc'][ieventhh+padh][0])
                    bx.scatter(h['ycmc'][ieventhh+padh][0], h['xcmc'][ieventhh+padh][0], s=50, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                    bx.arrow(h['ycmc'][ieventhh+padh][0],h['xcmc'][ieventhh+padh][0],dymc,dxmc,head_width=5, facecolor='blue', label="Real direction", alpha=0.5)
                
                bx.scatter(h['yc'][ieventhh+padh][0], h['xc'][ieventhh+padh][0], s=50, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                # limits = ax.axis()
                bx.arrow(h['yc'][ieventhh+padh][0],h['xc'][ieventhh+padh][0],dy,dx,head_width=5, facecolor='red', label="Rec direction", alpha=0.5)
                bx.add_patch(rect12);bx.add_patch(rect22);bx.add_patch(rect32)
                bx.text(1,0,tith,ha='left',va="bottom", wrap=True,family='serif', style='italic',bbox={'facecolor': 'white','edgecolor':'none','alpha': 0.3,'pad': 0.1,'boxstyle':'round4'}, transform=ax.transAxes)
                # if mc:
                #     bx.set_ylim([min([h['xcmc'][ieventgg+padg][0],h['xcmc'][ieventgg+padg][0]+dxmc,h['xc'][ieventhh+padh][0],h['xc'][ieventhh+padh][0]+dx,-130]), max([h['xcmc'][ieventgg+padg][0],h['xcmc'][ieventgg+padg][0]+dxmc,h['xc'][ieventhh+padh][0],h['xc'][ieventhh+padh][0]+dx,130])])
                #     bx.set_xlim([min(h['ycmc'][ieventgg+padg][0],h['ycmc'][ieventgg+padg][0]+dymc,h['yc'][ieventhh+padh][0],h['yc'][ieventhh+padh][0]+dy,-150), max(h['ycmc'][ieventgg+padg][0],h['ycmc'][ieventgg+padg][0]+dymc,h['yc'][ieventhh+padh][0],h['yc'][ieventhh+padh][0]+dy,150)])
                # bx.set_ylim([min([h['xc'][ieventgg+padg][0],h['xc'][ieventgg+padg][0]+dx,-130]), max([h['xc'][ieventgg+padg][0],h['xc'][ieventgg+padg][0]+dx,55])])
                # bx.set_xlim([min(h['yc'][ieventgg+padg][0],h['yc'][ieventgg+padg][0]+dy,-150), max(h['yc'][ieventgg+padg][0],h['yc'][ieventgg+padg][0]+dy,150)])

            plt.subplots_adjust(wspace=0, hspace=0)
            # plt.tight_layout()
            if mc:
                ax.set_xlabel("Gamma")
                ax.set_title(("nhit: %d priE: %.2f" %(g['nhit'][ieventgg+padg][0], g['priE'][ieventgg+padg][0])),fontsize=10)
            else:
                ax.set_title(("nhit: %d" %g['nhit'][ieventgg+padg][0]),fontsize=10)
            ax.axis('equal')
            # ax.autoscale(False)
            ax.invert_yaxis()


            if mc:
                bx.set_xlabel("Hadorn")
                bx.set_title(("nhit: %d priE: %.2f" %(h['nhit'][ieventhh+padh][0], h['priE'][ieventhh+padh][0])),fontsize=10)
            else:
                bx.set_title(("nhit: %d" % h['nhit'][ieventhh+padh][0]),fontsize=10)
            bx.axis('equal')
            # bx.autoscale(False)
            bx.invert_yaxis()
            plt.legend(loc="lower right", bbox_to_anchor=(1.63,0.85),framealpha=0.3)
            plt.subplots_adjust(bottom=0.2, left=0.1, right=0.9, top=0.8)
            fig.colorbar(a, ax=[ax, bx], shrink=0.7,label=('%s'%bbj)) #,orientation='horizontal', pad=0.18
        #g
        elif ifgd:
            if mode.get() == "image":
                vmin = cg[bbj][ieventgg].min()+abs(cg[bbj][ieventgg].min()/5)
                vmax = cg[bbj][ieventgg].max()
                if bbj == "npe" or bbj == "npe_c":
                    norm = colors.LogNorm(clip=True,vmin=vmax/500,vmax=vmax)
                else:
                    norm = colors.Normalize(vmin=vmin, vmax=vmax)

                plt.xticks(np.linspace(0,60,5),np.linspace(int(detx.min()),int(detx.max()),5))
                plt.yticks(np.linspace(0,52,5),np.linspace(int(dety.max()),int(dety.min()),5))

                #g
                a = plt.imshow(ndi.gaussian_filter(cg[bbj][ieventgg], sigma=float(sigma)),aspect='equal',cmap=cm,interpolation='none',norm=norm)

                xc,yc = rxy2ij(g['xc'][ieventgg+padg][0],g['yc'][ieventgg+padg][0])
                dx,dy = direction1(g['phi'][ieventgg+padg][0],g['theta'][ieventgg+padg][0])
                dx*-1;dy*=-1
                if mc:
                    xcmc,ycmc = rxy2ij(g['xcmc'][ieventgg+padg][0],g['ycmc'][ieventgg+padg][0])
                    dxmc,dymc = direction1(g['azimc'][ieventgg+padg][0],g['zenmc'][ieventgg+padg][0])
                    dymc*=-1
                    plt.scatter(xcmc, ycmc, s=50, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                    plt.arrow(xcmc,ycmc,dxmc,dymc,head_width=1, facecolor='blue', label="Real direction", alpha=0.5)
                plt.scatter(xc, yc, s=50, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                plt.arrow(xc,yc,dx,dy,head_width=1, facecolor='red', label="Rec direction", alpha=0.5)
                rect1 = plt.Rectangle((-0.5,-0.5),60,22,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect2 = plt.Rectangle((-0.5,22),30,29.5,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect3 = plt.Rectangle((29.5,22),30,29.5,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3)
                plt.gca().add_patch(rect1);plt.gca().add_patch(rect2);plt.gca().add_patch(rect3)
                plt.text(0,0,titg,ha='left',va="bottom",wrap=True,family='serif', style='italic',bbox={'facecolor': 'white','edgecolor':'none','alpha': 0.3,'pad': 0.1,'boxstyle':'round4'}, transform=plt.gca().transAxes)
                # if mc:
                #     plt.set_xlim([min([xcmc,xcmc+dxmc,xc,-5]), max([xcmc,xcmc+dxmc,xc,55])])
                #     plt.set_ylim([min(ycmc,ycmc+dymc,yc,yc+dy,-5), max(ycmc,ycmc+dymc,yc,yc+dy,65)])
                # plt.set_xlim([min([xc,xc+dx,-5]), max([xc,xc+dx,55])])
                # plt.set_ylim([min(yc,yc+dy,-5), max(yc,yc+dy,65)])

            else:
                vmin = g[bbj][ieventgg+padg][0].min()+abs(cg[bbj][ieventgg].min()/5)
                vmax = g[bbj][ieventgg+padg][0].max()
                if bbj == "npe" or bbj == "npe_c":
                    norm = colors.LogNorm(clip=True,vmin=vmax/500,vmax=vmax)
                else:
                    norm = colors.Normalize(vmin=vmin, vmax=vmax)
                a = plt.scatter(g['x'][ieventgg+padg][0], g['y'][ieventgg+padg][0], s=50, c=g[bbj][ieventgg+padg][0], alpha=0.5, norm=norm, cmap=cm)

                dx,dy = direction2(g['phi'][ieventgg+padg][0],g['theta'][ieventgg+padg][0])

                if mc:
                    dxmc,dymc = direction2(g['azimc'][ieventgg+padg][0],g['zenmc'][ieventgg+padg][0])
                    plt.scatter(g['xcmc'][ieventgg+padg][0], g['ycmc'][ieventgg+padg][0], s=60, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                    plt.arrow(g['xcmc'][ieventgg+padg][0],g['ycmc'][ieventgg+padg][0],dxmc,dymc,head_width=5, facecolor='blue', label="Real direction", alpha=0.5)
                plt.scatter(g['xc'][ieventgg+padg][0], g['yc'][ieventgg+padg][0], s=60, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                plt.arrow(g['xc'][ieventgg+padg][0],g['yc'][ieventgg+padg][0],dx,dy,head_width=5, facecolor='red', label="Rec direction", alpha=0.5)
                rect1 = plt.Rectangle((-150,-135),150,150,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect2 = plt.Rectangle((5,-135),150,150,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect3 = plt.Rectangle((-147.5,20),300,110,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3)
                plt.gca().add_patch(rect1);plt.gca().add_patch(rect2);plt.gca().add_patch(rect3)
                plt.text(0,0,titg,ha='left',va="bottom",wrap=True,family='serif', style='italic',bbox={'facecolor': 'white','edgecolor':'none','alpha': 0.3,'pad': 0.1,'boxstyle':'round4'}, transform=plt.gca().transAxes)
                # if mc:
                #     plt.set_xlim([min([xcmc,xcmc+dxmc,xc,-130]), max([xcmc,xcmc+dxmc,xc,130])])
                #     plt.set_ylim([min(ycmc,ycmc+dymc,yc,yc+dy,-150), max(ycmc,ycmc+dymc,yc,yc+dy,150)])
                # plt.set_xlim([min([xc,xc+dx,-130]), max([xc,xc+dx,55])])
                # plt.set_ylim([min(yc,yc+dy,-150), max(yc,yc+dy,150)])
            plt.axis('equal')
            if mc:
                plt.xlabel("Gamma")
                plt.title(("nhit: %d priE: %.2f" %(g['nhit'][ieventgg+padg][0], g['priE'][ieventgg+padg][0])),fontsize=10)
            else:
                plt.title(("nhit: %d" %g['nhit'][ieventgg+padg][0]),fontsize=10)
            plt.legend(loc="lower right", bbox_to_anchor=(1.33,0.85),framealpha=0.3)
            fig.colorbar(a, shrink=0.7,label=('%s'%bbj)) #, orientation='horizontal', pad=0.18
        #h
        elif ifhd:
            if mode.get() == "image":
                vmin = ch[bbj][ieventhh].min()+abs(ch[bbj][ieventhh].min()/5)
                vmax = ch[bbj][ieventhh].max()
                if bbj == "npe" or bbj == "npe_c":
                    norm = colors.LogNorm(clip=True,vmin=vmax/500,vmax=vmax)
                else:
                    norm = colors.Normalize(vmin=vmin, vmax=vmax)
                plt.xticks(np.linspace(0,60,5),np.linspace(int(detx.min()),int(detx.max()),5))
                plt.yticks(np.linspace(0,52,5),np.linspace(int(dety.max()),int(dety.min()),5))
                b = plt.imshow(ndi.gaussian_filter(ch[bbj][ieventhh], sigma=float(sigma)),aspect='equal',cmap=cm,interpolation='none',norm=norm)

                xc,yc = rxy2ij(h['xc'][ieventhh+padh][0],h['yc'][ieventhh+padh][0])
                dx,dy = direction1(h['phi'][ieventhh+padh][0],h['theta'][ieventhh+padh][0])
                dx*-1;dy*=-1

                if mc:
                    xcmc,ycmc = rxy2ij(h['xcmc'][ieventhh+padh][0],h['ycmc'][ieventhh+padh][0])
                    dxmc,dymc = direction1(h['azimc'][ieventhh+padh][0],h['zenmc'][ieventhh+padh][0])
                    dymc*=-1
                    plt.scatter(xcmc, ycmc, s=50, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core",marker="P")
                    plt.arrow(xcmc,ycmc,dxmc,dymc,head_width=1, facecolor='blue', label="Real direction", alpha=0.5)

                plt.scatter(xc, yc, s=50, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core",marker="X")
                plt.arrow(xc,yc,dx,dy,head_width=1, facecolor='red', label="Rec direction", alpha=0.5)
                rect1 = plt.Rectangle((-0.5,-0.5),60,22,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect2 = plt.Rectangle((-0.5,22),30,29.5,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect3 = plt.Rectangle((29.5,22),30,29.5,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3)
                plt.gca().add_patch(rect1);plt.gca().add_patch(rect2);plt.gca().add_patch(rect3)
                plt.text(0,0,tith,ha='left',va="bottom",wrap=True,family='serif', style='italic',bbox={'facecolor': 'white','edgecolor':'none','alpha': 0.3,'pad': 0.1,'boxstyle':'round4'}, transform=plt.gca().transAxes)
                # if mc:
                #     plt.set_xlim([min([xcmc,xcmc+dxmc,xc,-5]), max([xcmc,xcmc+dxmc,xc,55])])
                #     plt.set_ylim([min(ycmc,ycmc+dymc,yc,yc+dy,-5), max(ycmc,ycmc+dymc,yc,yc+dy,65)])
                # plt.set_xlim([min([xc,xc+dx,-5]), max([xc,xc+dx,55])])
                # plt.set_ylim([min(yc,yc+dy,-5), max(yc,yc+dy,65)])
            else:
                vmin = h[bbj][ieventhh+padh][0].min()+abs(ch[bbj][ieventhh].min()/5)
                vmax = h[bbj][ieventhh+padh][0].max()
                if bbj == "npe" or bbj == "npe_c":
                    norm = colors.LogNorm(clip=True,vmin=vmax/500,vmax=vmax)
                else:
                    norm = colors.Normalize(vmin=vmin, vmax=vmax)
                b = plt.scatter(h['x'][ieventhh+padh][0], h['y'][ieventhh+padh][0], s=50, c=h[bbj][ieventhh+padh][0], alpha=0.5, norm=norm, cmap=cm)
                dx,dy = direction2(h['phi'][ieventhh+padh][0],h['theta'][ieventhh+padh][0])

                if mc:
                    dxmc,dymc = direction2(h['azimc'][ieventhh+padh][0],h['zenmc'][ieventhh+padh][0])
                    plt.scatter(h['xcmc'][ieventhh+padh][0], h['ycmc'][ieventhh+padh][0], s=50, c='blue', alpha=0.8, norm=norm, cmap=cm,label="Real core:",marker="P")
                    plt.arrow(h['xcmc'][ieventhh+padh][0],h['ycmc'][ieventhh+padh][0],dxmc,dymc,head_width=5, facecolor='blue', label="Real direction", alpha=0.5)

                plt.scatter(h['xc'][ieventhh+padh][0], h['yc'][ieventhh+padh][0], s=50, c='red', alpha=0.8, norm=norm, cmap=cm, label="Rec core:",marker="X")
                plt.arrow(h['xc'][ieventhh+padh][0],h['yc'][ieventhh+padh][0],dx,dy,head_width=5, facecolor='red', label="Rec direction", alpha=0.5)
                rect1 = plt.Rectangle((-150,-135),150,150,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect2 = plt.Rectangle((5,-135),150,150,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3);rect3 = plt.Rectangle((-147.5,20),300,110,alpha=1,facecolor='none',edgecolor='black',linestyle='-.',linewidth=0.3)
                plt.gca().add_patch(rect1);plt.gca().add_patch(rect2);plt.gca().add_patch(rect3)
                plt.text(0,0,tith,ha='left',va="bottom",wrap=True,family='serif', style='italic',bbox={'facecolor': 'white','edgecolor':'none','alpha': 0.3,'pad': 0.1,'boxstyle':'round4'}, transform=plt.gca().transAxes)
                # if mc:
                #     plt.set_xlim([min([xcmc,xcmc+dxmc,xc,-130]), max([xcmc,xcmc+dxmc,xc,130])])
                #     plt.set_ylim([min(ycmc,ycmc+dymc,yc,yc+dy,-150), max(ycmc,ycmc+dymc,yc,yc+dy,150)])
                # plt.set_xlim([min([xc,xc+dx,-130]), max([xc,xc+dx,55])])
                # plt.set_ylim([min(yc,yc+dy,-150), max(yc,yc+dy,150)])
            if mc:
                plt.xlabel("Hadorn")
                plt.title(("nhit: %d priE: %.2f" %(h['nhit'][ieventhh+padh][0], h['priE'][ieventhh+padh][0])),fontsize=10)
            else:
                plt.title(("nhit: %d" % h['nhit'][ieventhh+padh][0]),fontsize=10)
            plt.axis('equal')
            plt.legend(loc="lower right", bbox_to_anchor=(1.33,0.85),framealpha=0.3)
            fig.colorbar(b, shrink=0.7,label=('%s'%bbj))
        else:
            BTt.insert(INSERT,"Error: \n");BTt.see(END)
            BTt.insert(INSERT,"Please choose a file.\n");BTt.see(END)

        titg = ""
        tith = ""
        for item in Parameterso:
            if not mc:
                if "mc" in item or item=='priE':
                    continue
            if CV[item].get():
                if item == "pinc":
                    if ifgd:
                        titg += str(item)+(("2:%.2f "%g[str(item)][ieventgg+padg][0][2]) + "\t")
                    if ifhd:
                        tith += str(item)+(("2:%.2f "%h[str(item)][ieventhh+padh][0][2]) +"\t")
                elif item == "nhit" or item == "evt":
                    if ifgd:
                        titg += str(item)+((":%d "%g[str(item)][ieventgg+padg][0]) + "\t")
                    if ifhd:
                        tith += str(item)+((":%d " %h[str(item)][ieventhh+padh][0]) + "\t")
                elif item == "mjd":
                    if ifgd:
                        tg = t0+datetime.timedelta(days=g[str(item)][ieventgg+padg][0]+1)
                        tg = datetime.datetime.date(tg)
                        titg += str(item)+((":%s " %tg) +"\n")
                    if ifhd:
                        th = t0+datetime.timedelta(days=h[str(item)][ieventgg+padg][0]+1)
                        th = datetime.datetime.date(th)
                        tith += str(item)+((":%s " %th) +"\n")
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
        canvas_spice.draw()
        toolbar.update()
        # plt.show()
        # if ss:
        #     print("c")
        #     canvas_spice2.draw()
        #     toolbar2.update()
        #     newWindow.update()

        # plt.ion()
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

def plotw():
    # global ss
    # global frame222
    # global frame2222    
    # global newWindow
    # global root
    # global fig
    # fig = plt.figure(figsize=(6.6,3.95),dpi=80) 
    # ss=1
    # newWindow = Toplevel(root,height=450, width=500)
    # frame222 = Frame(newWindow)
    # frame222.grid()
    # canvas_spice2 = FigureCanvasTkAgg(fig,frame222)
    # frame2222 = Frame(frame222)
    # frame2222.grid(row=1)
    # toolbar2 = NavigationToolbar2Tk(canvas_spice2, frame2222)
    # canvas_spice2.get_tk_widget().grid(row=0,column=0,columnspan=2)
    # canvas_spice2._tkcanvas.grid(row=0,column=0,columnspan=2)      #.pack(side=TOP, fill=BOTH, expand=1)
    # Button(frame222, text='Next>', command=plot).grid(row=1, column=1)
    # canvas_spice2.draw()
    # toolbar2.update()
    # root.update()
    plt.show()
    # plt.waitforbuttonpress()
    return 0

def load():
    global l1
    global p1
    global pltobd
    global Parameters
    global Parameterso
    global Parametersoo
    global ifgg
    global ifhh
    global fgname
    global fhname
    global tgamma
    global thadron
    global dg
    global dh
    global mc

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
    global mc
    try:
        g['%s'%item] = tree2array(tgamma, branches=['%s' %item])
        p1['value']+=30
        root.update()
        BTt.insert(INSERT,"Load "+item+" of gamma: Done!!!"+"\n");BTt.see(END)
        if ("mc" in item) or ( item == "priE"):
            mc=1
    except Exception as e:
        l1['text']=('no %s' %item)
        BTt.insert(INSERT,"ERROR:"+"\n");BTt.see(END)
        BTt.insert(INSERT,str(e.args)+"\n");BTt.see(END)
        BTt.insert(INSERT,"============="+"\n");BTt.see(END)
        BTt.insert(INSERT,traceback.format_exc()+"\n");BTt.see(END)
        if ("mc" in item) or ( item == "priE"):
            mc=0

def hload(item):
    global mc
    try:
        h['%s'%item] = tree2array(thadron, branches=['%s' %item])
        p1['value']+=30
        root.update()
        BTt.insert(INSERT,"Load "+item+" of hadron: Done!!!"+"\n");BTt.see(END)
        if ("mc" in item) or ( item == "priE"):
            mc=1
    except Exception as e:
        l1['text']=('no %s' %item)
        BTt.insert(INSERT,"ERROR:"+"\n");BTt.see(END)
        BTt.insert(INSERT,str(e.args)+"\n");BTt.see(END)
        BTt.insert(INSERT,"============="+"\n");BTt.see(END)
        BTt.insert(INSERT,traceback.format_exc()+"\n");BTt.see(END)
        if ("mc" in item) or ( item == "priE"):
            mc=0

def dgload(item):
    global mc
    try:
        g['%s'%item] = dg['%s'%item]
        p1['value']+=30
        BTt.insert(INSERT,"Load "+item+" of gamma: Done!!!"+"\n");BTt.see(END)
        root.update()
        if ("mc" in item) or ( item == "priE"):
            mc=1
    except Exception as e:
        l1['text']=('no %s' %item)
        BTt.insert(INSERT,"ERROR:"+"\n");BTt.see(END)
        BTt.insert(INSERT,str(e.args)+"\n");BTt.see(END)
        BTt.insert(INSERT,"============="+"\n");BTt.see(END)
        BTt.insert(INSERT,traceback.format_exc()+"\n");BTt.see(END)
        if ("mc" in item) or ( item == "priE"):
            mc=0

def dhload(item):
    global mc
    try:
        h['%s'%item] = dh['%s'%item]
        p1['value']+=30
        BTt.insert(INSERT,"Load "+item+" of hadron: Done!!!"+"\n");BTt.see(END)
        root.update()
        if ("mc" in item) or ( item == "priE"):
            mc=1
    except Exception as e:
        l1['text']=('no %s' %item)
        BTt.insert(INSERT,"ERROR:"+"\n");BTt.see(END)
        BTt.insert(INSERT,str(e.args)+"\n");BTt.see(END)
        BTt.insert(INSERT,"============="+"\n");BTt.see(END)
        BTt.insert(INSERT,traceback.format_exc()+"\n");BTt.see(END)
        if ("mc" in item) or ( item == "priE"):
            mc=0
    
if __name__ == "__main__":
    # Parameters = ['x', 'y', 'npe', 't', 'tr', 'dr', 'cpx', 'pinc', 'priE', 'nhit']
    Parameters = ['x','y','npe','t','xcmc','ycmc','irun','evt','nhit','nfit','sump1','sump2','mjd','sigma','theta','phi','xc','yc','ra','dec','fee','ch','iconf','npe_c','dr','tr','cpx','pinc','csr','dcore','dangle','priE','zenmc','azimc','h1','weight','dangle_real','dcore_real']
    Parametersl = ['x', 'y', 'fee', 'ch', 'iconf', 'npe', 'npe_c', 'dr', 'tr', 't']
    Parametersll = ['npe', 'npe_c', 'tr', 't', 'dr', 'fee', 'ch', 'iconf']
    Parametersm = ['x','y','xcmc','ycmc','xc','yc','zenmc','azimc','theta','phi',"nhit",'priE','evt']
    Parametersu = ['npe', 'npe_c', 't', 'tr', 'dr', 'fee', 'ch', 'iconf', 'irun', 'nfit', 'sump1', 'sump2', 'sigma','cpx', 'pinc', 'mjd', 'dcore', 'dangle', 'dcore_real', 'dangle_real', 'ra', 'dec', 'csr', 'h1', 'weight']
    # Parameterso = ['theta','phi','zenmc','azimc','irun', 'evt', 'nfit', 'sump1', 'sump2', 'cpx', 'pinc',  'mjd', 'dcore', 'dangle','dcore_real', 'dangle_real', 'ra', 'dec', 'sigma', 'csr', 'h1', 'weight']
    Parameterso = ['xc', 'yc', 'xcmc', 'ycmc', 'priE', 'nhit', 'theta','phi','zenmc','azimc', 'cpx', 'pinc', 'irun', 'evt', 'nfit', 'sump1', 'sump2', 'mjd', 'dcore', 'dangle','dcore_real', 'dangle_real', 'ra', 'dec', 'sigma', 'csr', 'h1', 'weight']
    Parametersoo = ['evt', 'cpx', 'pinc', 'theta','zenmc','nfit', 'sump1', 'sump2', 'mjd', 'dcore', 'dangle','dcore_real', 'dangle_real', 'ra', 'dec', 'sigma', 'csr', 'h1', 'irun', 'weight']
    # pl = len(list(set(Parameters) ^ set(Parametersm)))
    pl = len(Parametersu)
    global scc
    global mc
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
    ifrg = 0
    ifrh = 0
    mc = 0

    root = Tk()
    ww=500
    hh=450
    root.title('WCDA-Viewer V1.4.20')
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

    Label(frame1, text='The ROOT/npz file of Gamma/file1:   ').grid(row=1, column=0, sticky=W)

    e1 = Combobox(frame1, values=['/eos/user/z/zwang/wcda/data/mc/g4s2_rec/gamma_rec_143800t.root', 'datag.npz', './ES.146229.WCDA_EVENT.P111MC30TH5400_M100_Z.es-2.20220608071210.2806.reduced_pool5_rec.root'], width=37)
    e1.grid(row=1, column=1)
    Label(frame1, text='The ROOT/npz file of Hadron/file2:  ').grid(row=2, column=0, sticky=W)
    e2 = Combobox(frame1, values=['/eos/user/z/zwang/wcda/data/mc/g4s2_rec/hadron_rec_1970s.root', 'datah.npz', './ES.146229.WCDA_EVENT.P111MC30TH5400_M100_Z.es-2.20220608071210.2806.reduced_pool5_rec.root'], width=37)
    e2.grid(row=2, column=1)

    global C
    global CV
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
    bg = Checkbutton(frame4, text = "Gamma/file1?", onvalue = 1, offvalue = 0, width=13, variable=ifg)
    bh = Checkbutton(frame4, text = "Hadron/file2?", onvalue = 1, offvalue = 0, width=13, variable=ifh)
    bg.grid(row=3, column=1)
    bh.grid(row=3, column=2)


    global sc
    sc = StringVar(value="( nhit >= 1000 )")
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
    global fig2
    fig = plt.figure(figsize=(6.6,3.95),dpi=80)
    plt.close()

    canvas_spice = FigureCanvasTkAgg(fig,frame2)
    frame22 = Frame(frame2)
    frame22.grid(row=1)
    toolbar = NavigationToolbar2Tk(canvas_spice, frame22)
    canvas_spice.get_tk_widget().grid(row=0,column=0,columnspan=3)
    canvas_spice._tkcanvas.grid(row=0,column=0,columnspan=3)      #.pack(side=TOP, fill=BOTH, expand=1)
    Button(frame2, text='^popup^', command=plotw).grid(row=1, column=1,sticky=W)
    Button(frame2, text='Next>', command=plot).grid(row=1, column=2,sticky=W)

    # global toolbar2
    # global canvas_spice2
    # global newWindow
    # newWindow=None
    # frame2222=None
    # canvas_spice2=FigureCanvasTkAgg(fig,newWindow)
    # toolbar2=NavigationToolbar2Tk(canvas_spice, frame2222)
    # global ss
    # ss=0



    print('Show windows')
    root.protocol("WM_DELETE_WINDOW", lambda: sys.exit(0))
    root.mainloop()