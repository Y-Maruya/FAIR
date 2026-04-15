#!/bin/env python3

import os,copy,sys,datetime
import ROOT
from collections import defaultdict
from multiprocessing import Process, Semaphore

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection,Line3DCollection
import matplotlib.cm as cm
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator
from matplotlib import colormaps 

matplotlib.use("Agg")
np.random.seed(42)

plt.rcParams.update({
    "figure.figsize": (6, 6),
    "figure.dpi": 100,
    "savefig.dpi": 100
})

adcMin = 400
adcMax = 1200
def getCmap():
    norm = colors.Normalize(vmin=adcMin, vmax=adcMax)
    cmap = colormaps["rainbow"]
    return lambda adc: cmap(norm(adc))

z_layer = 1600/40 # TODO: approx
class reconstruction:
    def __init__(self,data,adc_key="hg_adc"):

        # Thanks to Yasu
        self._Pos_X=[100.2411,100.2411,100.2411,59.94146,59.94146,59.94146,19.64182,19.64182,19.64182,19.64182,59.94146,100.2411,100.2411,59.94146,19.64182,100.2411,59.94146,19.64182,-20.65782,-60.95746,-101.2571,-20.65782,-60.95746,-101.2571,-101.2571,-60.95746,-20.65782,-20.65782,-20.65782,-20.65782,-60.95746,-60.95746,-60.95746,-101.2571,-101.2571,-101.2571]
        self._Pos_Y=[141.04874,181.34838,221.64802,141.04874,181.34838,221.64802,141.04874,181.34838,221.64802,261.94766,261.94766,261.94766,302.2473,302.2473,302.2473,342.54694,342.54694,342.54694,342.54694,342.54694,342.54694,302.2473,302.2473,302.2473,261.94766,261.94766,261.94766,221.64802,181.34838,141.04874,221.64802,181.34838,141.04874,221.64802,181.34838,141.04874]
        self.chip_dis_X=239.3
        self.chip_dis_Y=241.8
        self.HBU_X=239.3
        self.HBU_Y=725.4
        global z_layer
        self.z_layer = z_layer

        nHits = data["nHits"]
        self.hitPositions = []
        cmap = getCmap()
        for iHit in range(len(data["hittag"])):
            if not data["hittag"][iHit]: continue
            layer = data["layer"][iHit]
            channel = data["channel"][iHit]
            asic = data["asic"][iHit]
            adc = data[adc_key][iHit]

            color = cmap(adc)

            x = self.Pos_X(channel,asic)
            y = self.Pos_Y(channel,asic)
            z = self.Pos_Z(layer)
            d = {"x":x,"y":y,"z":z,"adc":adc,"color":color}
            self.hitPositions.append(d)


    def Pos_X(self,channel,asic):
        HBU_ID=asic/3;
        asic=asic%3;
        # if(asic!=0){
        #     if(channel==2)channel=0;
        #     else if(channel==0)channel=2;
        #     if(channel==33)channel=35;
        #     else if(channel==35)channel=33;
        # }
        return 500+(self._Pos_Y[channel]-asic*self.chip_dis_Y);


    def Pos_Y(self,channel,asic):
        HBU_ID=asic/3;
        return 500-(-self._Pos_X[channel]+(HBU_ID-1)*self.HBU_X);

    def Pos_Z(self,layer):
        return layer*self.z_layer

class eventdisplay:
    def __init__(self,iPath,oDir,run="",events=[],filters=defaultdict(bool),event_counter=None):
        self.iPath = iPath
        self.oDir = oDir
        self.run = run
        self.events = events
        self.filters = filters
        self.event_counter = event_counter
        self.loadData()
        self.sema = Semaphore(1)

        self.xExt = 1000
        self.yExt = 1000
        self.zExt = 1600

    def loadData(self):
        iFile = ROOT.TFile.Open(self.iPath)
        iTree = iFile.Get("events")
        nEvents = iTree.GetEntries()
        self.data = defaultdict(dict)
        print(f"Total events:{nEvents}")

        if self.event_counter is not None:
            target_entry = None
            for event in range(nEvents):
                iTree.GetEntry(int(event))
                if int(getattr(iTree,"EventSummary.event_counter")) == int(self.event_counter):
                    target_entry = event
                    break
            if target_entry is None:
                print(f"Event with event_counter {self.event_counter} not found")
                iFile.Close()
                self.events = []
                return
            self.events = [target_entry]
        else:
            self.events = set(np.random.randint(0,nEvents-1,self.events))

        for event in self.events:
            if event>=nEvents:
                print(f"Skipping event {event}")
                continue
            iTree.GetEntry(int(event))
            nHits = getattr(iTree,"EventSummary.nHits")
            self.data[event]["nHits"] = nHits
            # Get trigger inputs from TLURawData.Inputs vector
            triggerInputs = list(getattr(iTree,"TLURawData.Inputs"))
            print(f"Event {event}: nHits={nHits}, triggerInputs={triggerInputs}")
            # # filtering events
            # # -- select events with 2 scintillator hits
            # if self.filters["veto2"] and (len(triggerInputs) < 6 or not triggerInputs[4] or not triggerInputs[5]): 
            #     continue
            # # -- select events with 0 scintillator hits
            # if self.filters["veto0"] and (len(triggerInputs) >= 4 and (triggerInputs[4] or triggerInputs[5])): 
            #     continue
            # # -- select events with nHits
            # if self.filters["nHit"] and nHits<self.filters["nHit"]: 
            #     continue

            # # Store event summary data
            # for attr in dir(getattr(iTree, "EventSummary")):
            #     if not attr.startswith("_"):
            #         try:
            #             self.data[event][f"EventSummary.{attr}"] = getattr(getattr(iTree, "EventSummary"), attr)
            #         except:
            #             pass
            
            # # Store TLURawData
            # for attr in dir(getattr(iTree, "TLURawData")):
            #     if not attr.startswith("_") and attr not in ["Inputs", "FineTimestamps"]:
            #         try:
            #             self.data[event][f"TLURawData.{attr}"] = getattr(getattr(iTree, "TLURawData"), attr)
            #         except:
            #             pass
            
            # Store vector data
            self.data[event]["TLURawData.Inputs"] = triggerInputs
            self.data[event]["TLURawData.FineTimestamps"] = list(getattr(iTree,"TLURawData.FineTimestamps"))
            
            # Store RawHits data (the actual detector hits)
            self.data[event]["RawHits.v.cellID"] = list(getattr(iTree, "RawHits.v.cellID"))
            self.data[event]["RawHits.v.hg_adc"] = list(getattr(iTree, "RawHits.v.hg_adc"))
            self.data[event]["RawHits.v.lg_adc"] = list(getattr(iTree, "RawHits.v.lg_adc"))
            self.data[event]["RawHits.v.hittag"] = list(getattr(iTree, "RawHits.v.hittag"))
            self.data[event]["RawHits.v.bcid"] = list(getattr(iTree, "RawHits.v.bcid"))

            # Store RecoHits data
            # recoHits = getattr(iTree, "RecoHits")
            # self.data[event]["RecoHits.v.cellID"] = list(recoHits.v.cellID)
            # self.data[event]["RecoHits.v.Edep"] = list(recoHits.v.Edep)
            # self.data[event]["RecoHits.v.Nmip"] = list(recoHits.v.Nmip)
            
        iFile.Close()
        self.events = list(self.data.keys())



    def draw(self):
        ps = []
        batches = []
        batchCounter = 0
        for event in self.events:
            batchCounter+=1
            ps.append(Process(target=self.drawEvent,args=(event,self.data[event],)))
            if batchCounter>100:
                batchCounter=0
                batches.append(ps)
                ps = []
        batches.append(ps)


        for iBatch,batch in enumerate(batches):
            print("Batch",iBatch)
            for p in batch: p.start()
            for p in batch: p.join()
            for p in batch: p.close()

    def drawHit(self,hit):
        x = hit["x"]
        y = hit["y"]
        z = hit["z"]
        adc = hit["adc"]
        color = hit["color"]
        dx = 50
        dy = 50
        verts = [
            [x+0,  y+0,  z],
            [x+dx, y+0,  z],
            [x+dx, y+dy, z],
            [x+0,  y+dy, z]
        ]
        rect = Poly3DCollection([verts],color=color)
        return rect

    def drawWall(self,x,y,z,dx,dy,rot=0,color="C0",outline=True,alpha=0.1):
        verts = np.array([
            [x+0,  y+0,  z],
            [x+dx, y+0,  z],
            [x+dx, y+dy, z],
            [x+0,  y+dy, z]
        ])
        for i in range(rot):
            verts = verts[:, [1, 2, 0]]

        ret = []
        rect = Poly3DCollection([verts],alpha=alpha,color=color)
        ret.append(rect)

        if outline:
            verts = np.vstack([verts, verts[0]])
            rect = Line3DCollection([verts],color="k",lw=0.5)
            ret.append(rect)

        return ret

    def drawAhcal(self,triggerInputs):

        xExt = self.xExt
        yExt = self.yExt
        zExt = self.zExt

        ret = []
        boxColor = "C0"
        boxAlpha = 0.3
        ret+=self.drawWall(0,0,0,xExt,yExt,color=boxColor,alpha=boxAlpha) # downstream
        ret+=self.drawWall(0,0,zExt,xExt,yExt,color=boxColor,alpha=boxAlpha) # upstream
        ret+=self.drawWall(0,0,0,zExt,yExt,rot=1,color=boxColor,alpha=boxAlpha) # top
        ret+=self.drawWall(0,0,xExt,zExt,yExt,rot=1,color=boxColor,alpha=boxAlpha) # bottom
        ret+=self.drawWall(0,0,0,xExt,zExt,rot=2,color=boxColor,alpha=boxAlpha) # back
        # ret+=self.drawWall(0,0,xExt,xExt,zExt,rot=2) # front

        # trigger layers
        self.triggerLayers = [9,19,29,38] # TODO: correct layers
        for iTrigger,triggerLayer in enumerate(self.triggerLayers):
            zPos = z_layer*triggerLayer
            color = "y" if triggerInputs[iTrigger] else "#eee"
            alpha = 0.15 if triggerInputs[iTrigger] else 0.15
            ret+=self.drawWall(0,0,zPos,xExt,yExt,color=color,outline=0,alpha=alpha)

        # scintillators
        dScint = 200
        for scint in [4,5]:
            color = "r" if triggerInputs[scint] else "#eee"
            alpha = 0.22 if triggerInputs[scint] else 0.25
            zPos = -dScint-(scint-4)*100
            ret+=self.drawWall(0,0,zPos,xExt,yExt,color=color,outline=0,alpha=alpha) 

        return ret

    def drawBackground(self,ax):

        nx, ny = 400, 400
        grad = np.linspace(0, 1, ny)
        grad = np.tile(grad, (nx, 1))

        ax.imshow(
            grad,
            # vmin = 0.5,
            # vmax = 1,
            vmin = -1,
            vmax = 0.5,
            extent=[-0.5,0.5,-0.5,0.5],
            origin='lower',
            cmap='Grays',
            aspect='auto',
            zorder=-100
        )

    def drawEvent(self,event,data):
        self.sema.acquire()
        print(f"Drawing event {event}")
        
        # Get trigger inputs from stored data
        triggerInputs = data.get("TLURawData.Inputs", [False]*6)
        # Ensure we have 6 inputs (veto channels are at indices 4 and 5)
        while len(triggerInputs) < 6:
            triggerInputs.append(False)

        # Prepare hit data for reconstruction
        hitData = {
            "layer": [],
            "channel": [],
            "asic": [],
            "hittag": data.get("RawHits.v.hittag", []),
            "hg_adc": data.get("RawHits.v.hg_adc", []),
            "lg_adc": data.get("RawHits.v.lg_adc", []),
            "nHits": data.get("nHits", 0)
        }
        
        # Extract layer, channel, asic from cellID
        # cellID encoding: layer*100000 + asic*10000 + channel
        for cellID in data.get("RawHits.v.cellID", []):
            layer = cellID // 100000
            asic = (cellID // 10000) % 10
            channel = cellID % 10000
            hitData["layer"].append(layer)
            hitData["channel"].append(channel)
            hitData["asic"].append(asic)

        # Get timestamp from TLURawData
        timestamp = data.get("TLURawData.Timestamp", 0) / 1000000
        dt = datetime.datetime.fromtimestamp(timestamp)
        timestamp_str = dt.strftime("%d/%m/%y %H:%M:%S")

        for adc_key, adc_tag in [("hg_adc", "hg"), ("lg_adc", "lg")]:
            reco = reconstruction(hitData, adc_key=adc_key)

            ax = plt.figure().add_subplot(projection='3d')
            ax.dist = 15

            [ax.add_collection3d(a) for a in self.drawAhcal(triggerInputs)]

            for hit in reco.hitPositions:
                rect = self.drawHit(hit)
                if rect:
                    ax.add_collection3d(rect)

            self.fontColor = "#eee"
            y0= self.yExt*1.45
            m = self.yExt*0.11
            s = 9
            ax.text(0,y0,self.zExt,"FASER AHCAL",size=11,color=self.fontColor)
            ax.text(0,y0-m,self.zExt,f"Run {self.run}, event {event} ({adc_tag})",size=s,color=self.fontColor)
            # ax.text(0,y0-2*m,self.zExt,timestamp_str,size=s,color=self.fontColor)
            ax.text(0,y0-3*m,self.zExt,f"{len(reco.hitPositions)} hits",size=s,color=self.fontColor)

            ax.view_init(elev=30, azim=10,roll=93)
            ax.set_box_aspect((1, 1, 1))

            # ax.dist = 20
            # ax.set_proj_type("persp")

            ax.set_box_aspect((1, 1, 1))

            x_middle = np.mean([0,self.xExt])
            y_middle = np.mean([0,self.yExt])-100
            z_middle = np.mean([-450,self.zExt])
            radius = 0.47 * max([self.xExt,self.yExt,self.zExt])
            ax.set_xlim3d([x_middle - radius, x_middle + radius])
            ax.set_ylim3d([y_middle - radius, y_middle + radius])
            ax.set_zlim3d([z_middle - radius, z_middle + radius])

            ax.set_axis_off()

            self.drawBackground(ax)

            smallHeight = 15

            # inset 2D scatter x,y
            s = f"{smallHeight}%"
            ax2 = inset_axes(ax, width=s, height=s, loc="lower left", borderpad=1,
                             bbox_to_anchor=(0.05, 0.02, 1, 1),
                             bbox_transform=ax.transAxes,
                            )
            for hit in reco.hitPositions:
                ax2.plot(hit["x"], hit["y"],color=hit["color"],marker='s',ms=2)
            ax2.set_ylim([0,self.yExt])
            ax2.set_xlim([0,self.xExt])
            self.prepSmallPlot(ax2,title="Front",x="x",y="y")

            # inset 2D scatter y,z
            sx = f"{smallHeight*self.zExt/self.yExt}%"
            sy = f"{smallHeight}%"
            ax3 = inset_axes(ax, width=sx, height=sy, loc="lower center", borderpad=1,
                             bbox_to_anchor=(-0.04, 0.02, 1, 1), bbox_transform=ax.transAxes,
                            )
            for hit in reco.hitPositions:
                ax3.plot(-hit["z"], hit["y"],color=hit["color"],marker='s',ms=2)
            ax3.set_ylim([0,self.yExt])
            ax3.set_xlim([-self.zExt,0])
            self.prepSmallPlot(ax3,title="Side",x="z",y="y")

            # inset 2D scatter x,z
            sx = f"{smallHeight*self.zExt/self.yExt}%"
            sy = f"{smallHeight}%"
            ax4 = inset_axes(ax, width=sx, height=sy, loc="lower right", borderpad=1,
                             bbox_to_anchor=(-0.05, 0.02, 1, 1), bbox_transform=ax.transAxes,
                            )
            for hit in reco.hitPositions:
                ax4.plot(-hit["z"], hit["x"],color=hit["color"],marker='s',ms=2)
            ax4.set_ylim([0,self.xExt])
            ax4.set_xlim([-self.zExt,0])
            self.prepSmallPlot(ax4,title="Top",x="z",y="x")

            # ADC colorscale
            adcHeight = 2
            adcWidth = 14
            sy = f"{adcHeight}%"
            sx = f"{adcWidth}%"
            ax5 = inset_axes(ax, width=sx, height=sy, loc="lower right", borderpad=1,
                             bbox_to_anchor=(-0.05, 0.23, 1, 1), bbox_transform=ax.transAxes,
                            )
            data_cbar = np.linspace(adcMin, adcMax, 50).reshape(1, 50)
            plt.imshow(
                data_cbar,
                aspect="auto",
                cmap=colormaps["rainbow"],
                alpha=0.75,
                norm=colors.Normalize(vmin=adcMin, vmax=adcMax)
            )
            ax5.set_axis_off()
            self.prepSmallPlot(ax5,title="",x="",y="")
            ax5.set_title(f"{adc_tag.upper()} ADC value",fontsize=7,pad=-1,color=self.fontColor,loc="right")

            label = "-".join([f"{k}" for k,v in self.filters.items() if v])
            plt.savefig(f"{self.oDir}/aed-{self.run}-{event}-{label}-{adc_tag}.png",bbox_inches="tight",dpi=400)
            plt.close()


        self.sema.release()

    def prepSmallPlot(self,ax,title="",x="",y=""):
        smallColor  = "#222"
        smallLabels = 6
        ax.set_title(title,fontsize=smallLabels,pad=-2,color=self.fontColor)
        ax.set_xlabel(x,color=self.fontColor)
        ax.set_ylabel(y,color=self.fontColor)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.xaxis.label.set_size(smallLabels)
        ax.yaxis.label.set_size(smallLabels)
        ax.xaxis.labelpad = -2
        ax.yaxis.labelpad = -2
        ax.xaxis.set_minor_locator(MultipleLocator(100))
        ax.yaxis.set_minor_locator(MultipleLocator(100))
        ax.minorticks_on()
        ax.set_facecolor(smallColor)
        self.setAxisColor(ax,self.fontColor)
        ticksInside(ax)

    def setAxisColor(self,ax,color):
        # Tick colors
        ax.tick_params(axis="x", which="major", colors=color)
        ax.tick_params(axis="y", which="major", colors=color)
        ax.tick_params(axis="x", which="minor", colors=color)
        ax.tick_params(axis="y", which="minor", colors=color)

        # Spine (axis line) colors
        ax.spines["bottom"].set_color(color)
        ax.spines["top"].set_color(color)
        ax.spines["left"].set_color(color)
        ax.spines["right"].set_color(color)

def ticksInside(ax,removeXLabel=False,removeYLabel=False):
    """ Make atlas style ticks """
    ax.tick_params(labeltop=False, labelright=False)
    plt.xlabel(ax.get_xlabel(), horizontalalignment='right', x=1.0)
    plt.ylabel(ax.get_ylabel(), horizontalalignment='right', y=1.0)
    ax.tick_params(axis='y',direction="in",labelleft=not removeYLabel,left=1,right=1,which='both')
    ax.tick_params(axis='x',direction="in",labelbottom=not removeXLabel,bottom=1, top=1,which='both')



if __name__=="__main__":
    iPath = sys.argv[1]
    oDir = sys.argv[2]
    run = iPath.find("run") #run*
    run = iPath[run+3:run+8] if run>=0 else "run0000"

    # Check if specific event_counter is provided
    if len(sys.argv) > 3:
        # Specific event_counter mode
        event_counter = int(sys.argv[3])
        events = 1
        filters = defaultdict(bool)
        ed = eventdisplay(iPath,oDir,run=run,events=events,filters=filters,event_counter=event_counter)
        if ed.events:
            event = ed.events[0]
            print(f"Drawing event_counter {event_counter} from tree entry {event}")
            ed.drawEvent(event_counter, ed.data[event])
    else:
        # Default mode: multiple event sets
        events = 30
        filters = defaultdict(bool)
        ed = eventdisplay(iPath,oDir,run=run,events=events,filters=filters)
        ed.draw()

        events = 5
        filters = defaultdict(bool)
        filters["veto0"] = True
        ed = eventdisplay(iPath,oDir,run=run,events=events,filters=filters)
        ed.draw()

        events = 5
        filters = defaultdict(bool)
        filters["veto2"] = True
        ed = eventdisplay(iPath,oDir,run=run,events=events,filters=filters)
        ed.draw()

        events = 5
        filters = defaultdict(bool)
        filters["nHit"] = 50
        ed = eventdisplay(iPath,oDir,run=run,events=events,filters=filters)
        ed.draw()
