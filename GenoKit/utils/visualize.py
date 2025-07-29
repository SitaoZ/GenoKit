import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle, FancyArrowPatch
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

class visual_transcript(object):
    '''
    Visualize the uORF/dORF in transcript location (relative position) 
    Parameters:
     - transcript_location: a 1D list postion of transcript
     - cds_location: a 1D list CDS postion of transcript
     - xorf_location: a 2D list position of uORFs/dORFs
    Return:
     - a figure object
    '''
    def __init__(self, output, transcript_id, transcript_location, cds_location, xorf_location, feature_type):
        # init
        self.output = output                  # output prefix of figure
        self.t_id = transcript_id             
        self.transcript = transcript_location # 1D list
        self.CDS = cds_location               # 1D list
        self.xORFs = xorf_location            # 2D list
        self.feature_type = feature_type      # uorf types: type 1 2 3
        # secondary 
        self.xt = self.transcript[0] # real position 
        self.xc = self.CDS[0]
        # figure global variable
        self.xticks = transcript_location
        self.yt = 10
        self.yc = 10
        self.height = 1
        self.scale = 0.5 
        self.tip = 1
        
    def draw(self):
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        # xtick 
        label_unit = (self.xticks[-1]-self.xticks[0])//10
        xtick = [i for i in range(self.xticks[0], self.xticks[-1],label_unit)]
        ax.scatter(xtick, len(xtick)*[self.yt+2.5],marker = '|',color='grey')
        for i,text in enumerate(xtick):
            plt.annotate(text, xy=(xtick[i],self.yt+3), ha='center',color='grey')
        # 1.transcript 
        ax.add_patch(Rectangle((self.xt, self.yt + 0.5*self.scale * self.height), \
                    self.transcript[1] - self.xt, \
                    self.scale * self.height, color='#DA70D6')
                    ) # xy, width, height
        plt.text(self.transcript[1]*1.01, \
                 self.yt + 0.5*self.scale * self.height, self.t_id, \
                 fontsize='small')
        # 2.CDS, CDS overlap with transcript 
        ax.add_patch(Rectangle((self.xc, self.yt), \
                    self.CDS[1] - self.xc, \
                    self.height, color='#FF8C00'))
        plt.text((self.CDS[1] + self.CDS[0])/2, \
                  self.yc + 1.5, \
                  "CDS", fontsize='small')
        # 3.UTR5
        if self.CDS[0] != self.transcript[0]:
            # 5 primer utr exist 
            plt.text((self.xt + self.xc)/2, \
                      self.yc + 1.5, \
                      "5' UTR", fontsize='small')
        # 4.UTR3
        if self.CDS[1] != self.transcript[1]:
            # 3 primer utr exist 
            plt.text((self.transcript[1] + self.CDS[1])/2, \
                      self.yc + 1.5, \
                      "3' UTR", fontsize='small')
        # 5.uORF/dORF
        for i,xorf in enumerate(self.xORFs):
            i = i + 1
            xu = xorf[0]
            # yu = self.yc - 2 * i * self.height
            yu = self.yt - 2 - i + self.scale * self.height
            ax.add_patch(Rectangle((xu, yu), \
                         xorf[1] - xu, 0.8 * self.height, color='#43CD80'))
            plt.text(xorf[1]*1.01, yu, self.t_id + '.'+self.feature_type, fontsize='small')
        plt.xlim(self.transcript[0] - 10, self.transcript[1] + 10)
        plt.ylim(-2 * len(self.xORFs),10 + 2 * len(self.xORFs))
        plt.axis('off')
        plt.savefig(f'{self.output}.pdf')
        plt.show()


#a = visual_transcript('out', 'aaa',[1,1000],\
#                        [100,800],\
#                        [[30,90],[60,440]], 'type1')
#a.draw()



class visual_transcript_genome(object):
    '''
    Visualize a transcript on genome, (absolute position)
    Parameters:
     - output: output file path, str
     - transcript_id: transcript id, str
     - exon_location: exon location blocks in genome, 2D list 
     - cds_location: CDS location blocks in genome, 2D list 
     - xorf_location: uORF/dORF location blocks in genome, 2D list 
    Return:
     - a figure object
    '''
    def __init__(self, strand, output, transcript_id, exon_location, cds_location, xorf_location, feature_type):
        # init
        self.strand = strand               # strand in genome
        self.output = output               # output prefix of a figure
        self.transcript_id = transcript_id
        self.exon_location = exon_location # 2D list 
        # print('exon_location=',exon_location)
        self.cds_location = cds_location   # 2D list 
        # print('cds_location=',cds_location)
        self.xorf_location = xorf_location # 3D list  
        # print("xorf_location=", xorf_location)
        self.feature_type = feature_type   # uorf types: type 1 2 3
        # figure global variable
        self.height = 1
        self.yt = 10                       # relative position
        self.yc = 10
        self.scale = 0.8
        self.tip = 1
        self.xticks = exon_location
        # intron from exon location 
        self.transcript_intron_location = self.transcript_intron_locate()
    
    def transcript_intron_locate(self):
        '''
        intron of transcript in genome coordinate, same to genome location 
        -1 
        '''
        intron_list = []
        intron_start = 0
        intron_end = 0
        # no intron 
        if len(self.exon_location) == 1:
            return None
        # intron exist 
        else:
            for i,exon in enumerate(self.exon_location):
                i = i + 1
                if i == 1 :
                    # first exon 
                    intron_start = exon[1]-1 # genome location should add 1 exon[1]+1; visualization ignore this
                    intron_list.append(intron_start)
                    continue 
                elif i == len(self.exon_location):
                    # last exon 
                    intron_end = exon[0]-1
                    intron_list.append(intron_end)
                else:
                    # inter exon
                    intron_list.append(exon[0]-1)
                    intron_list.append(exon[1]-1)
            array = np.array(intron_list)
            intron_location_2d = array.reshape(-1, 2)
            #if self.strand == '+':
            #    genome_intron_location_2d = intron_location_2d
            #else:
            #    genome_intron_location_2d = intron_location_2d[::-1]
            return intron_location_2d

    def cds_intron_location(self):
        '''
        intron of cds in genome coordinate, same to genome location 
        '''
        intron_list = []
        intron_start = 0
        intron_end = 0
        # no intron 
        if len(self.cds_location) == 1:
            return None
        # intron exist 
        else:
            for i,exon in enumerate(self.cds_location):
                i = i + 1
                if i == 1 :
                    # first exon 
                    intron_start = exon[1]-1 # genome location should add 1 exon[1]+1; visualization ignore this
                    intron_list.append(intron_start)
                    continue
                elif i == len(self.cds_location):
                    # last exon 
                    intron_end = exon[0]-1
                    intron_list.append(intron_end)
                else:
                    # inter exon
                    intron_list.append(exon[0]-1)
                    intron_list.append(exon[1]-1)
            array = np.array(intron_list)
            intron_location_2d = array.reshape(-1, 2)
            #if self.strand == '+':
            #    genome_intron_location_2d = intron_location_2d
            #else:
            #    genome_intron_location_2d = intron_location_2d[::-1]
            return intron_location_2d 
    
    def xorf_intron_location(self):
        '''
         uorf/dorf intron location  
        '''
        intron_list_all = []
        intron_start = 0
        intron_end = 0
        # 2D list 
        for xorf in self.xorf_location:
            # wotyout intron 
            # print('xorf:', xorf)
            if len(xorf) == 1:
                intron_list_all.append([])
            else:
                # with intron 
                if self.strand == '-':
                    xorf = xorf[::-1] # reverse list 
                intron_list = []
                for i,exon in enumerate(xorf):
                    i = i + 1
                    if i == 1 :
                        # first exon 
                        intron_start = exon[1]-1 # genome location should add 1 exon[1]+1; visualization ignore this
                        intron_list.append(intron_start)
                        continue
                    elif i == len(xorf):
                        # last exon 
                        intron_end = exon[0]-1
                        intron_list.append(intron_end)
                    else:
                        # inter exon
                        intron_list.append(exon[0]-1)
                        intron_list.append(exon[1]-1)
                array = np.array(intron_list)
                intron_location_2d = array.reshape(-1, 2)
                intron_list_all.append(intron_location_2d)    
                #if self.strand == '+':
                #    genome_intron_location_2d = intron_location_2d
                #else:
                #    genome_intron_location_2d = intron_location_2d[::-1]
        return intron_list_all
        
    def get_arrow_ab(self, exon):
        """get arrow posA and posB"""
        if self.strand == '+':
            pos_a = exon[0]
            pos_b = exon[1]
        else:
            pos_a = exon[1]
            pos_b = exon[0]
        return pos_a, pos_b
         
    def draw(self):
        #fig = plt.figure(figsize=(6, 8), dpi=300)
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        # x_ticks
        # ax.set_xticks([i for i in range(self.xticks[0][0], self.xticks[-1][1],50)])
        label_unit = (self.xticks[-1][1]-self.xticks[0][0])//10
        xtick = [i for i in range(self.xticks[0][0], self.xticks[-1][1],label_unit)]
        ax.scatter(xtick, len(xtick)*[self.yt+2],marker = '|',color='grey')
        for i,text in enumerate(xtick):
            if i % 2 == 0:
                plt.annotate(text, xy=(xtick[i],self.yt+2.5), ha='center',color='grey')
        #ax.plot([self.xticks[0][0], self.xticks[-1][1]], \
        #                [self.yt + 2, self.yt + 2], \
        #                solid_capstyle='butt', dash_capstyle = 'butt', linewidth=0.7,color='grey')
        # 1.transcript intron line for transcript
        for intron in self.transcript_intron_location:
            if len(intron) != 0 :
                ax.plot([intron[0], intron[1]], \
                        [self.yt + self.scale * self.height, \
                        self.yt + self.scale * self.height], \
                        solid_capstyle='butt', dash_capstyle = 'butt', linewidth=0.7,color='grey')
                # plt.axvline(x=intron[1])
        # 2.transcript exon rectangle 
        for exon in self.exon_location:
            # Rectangle
            #ax.add_patch(Rectangle((exon[0], self.yt), \
            #            exon[1]-exon[0]-1, self.height, color='#FF8C00')) # xy, width(exon length), height
             
            # FancyArrowPatch; optional method
            pos_a, pos_b = self.get_arrow_ab(exon)
            arrow = FancyArrowPatch(
               posA=(pos_a, self.yt + self.scale * self.height),
               posB=(pos_b, self.yt + self.scale * self.height), 
               fc='#FF8C00', ec='#FF8C00', shrinkA=0.0, shrinkB=0.0,
               arrowstyle='simple, head_width=5, head_length=5, tail_width=5', capstyle='butt')
            ax.add_artist(arrow)
            # plt.text((exon[0]+exon[1])/2, self.yt, 'exon', fontsize='small')
        # 3.transcript text 
        plt.text(self.exon_location[-1][1] * 1.01, \
                 self.yt + self.scale * self.height, \
                 self.transcript_id, fontsize=5)
        
        # cds intron line, different with transcript
        # 1.cds intron exist 
        if self.cds_intron_location().all():
            for intron in self.cds_intron_location():
                if len(intron) != 0 :
                    ax.plot([intron[0], intron[1]], \
                            [self.yt - 2 + self.scale * self.height, \
                             self.yt - 2 + self.scale * self.height], \
                             solid_capstyle='butt', dash_capstyle = 'butt', linewidth=0.7, color='grey')
        # 2.cds exon rectangle
        for cds in self.cds_location:
            # Rectangle
            #ax.add_patch(Rectangle((cds[0], self.yt - 2), \
            #            cds[1] - cds[0], self.scale * self.height, color='#DA70D6')) # xy, width(exon length), height  
            # FancyArrowPatch
            pos_a, pos_b = self.get_arrow_ab(cds)
            arrow = FancyArrowPatch(
                 posA=(pos_a, self.yt - 2 + self.scale * self.height), shrinkA=0.0, shrinkB=0.0,
                 posB=(pos_b, self.yt - 2 + self.scale * self.height), fc='#DA70D6', ec='#DA70D6',
                 arrowstyle='simple, head_width=5, head_length=5, tail_width=5', capstyle='butt')
            ax.add_artist(arrow)
        # 3.cds text 
        plt.text(self.cds_location[-1][1] * 1.01, \
                 self.yt - 2 + self.scale * self.height, \
                 self.transcript_id + '.CDS', fontsize=5)
        
        # uORF/dORF
        for i,xorf in enumerate(self.xorf_location):
            i = i + 1
            # 1.xorf intron line, different with transcript
            if len(self.xorf_intron_location()[i-1]) != 0:
                for intron in self.xorf_intron_location()[i-1]:
                    # intron 1D list 
                    if len(intron) != 0 :
                        ax.plot([intron[0], intron[1]], \
                            [self.yt - 2 - i + self.scale * self.height , \
                             self.yt - 2 - i  + self.scale * self.height], \
                            solid_capstyle='butt', dash_capstyle = 'butt', linewidth=0.7,color='grey')
                    else:
                        continue
            # 2.xrof exon rectangle
            for j,u in enumerate(xorf):
                pass
                # Rectangle
                #ax.add_patch(Rectangle((u[0], self.yt - 2 - i), \
                #           u[1]-u[0], self.scale * self.height,capstyle='butt', color='#43CD80'))
                pos_a, pos_b = self.get_arrow_ab(u)
                arrow = FancyArrowPatch(
                     posA=(pos_a, self.yt - 2 - i + self.scale * self.height), shrinkA=0.0, shrinkB=0.0,
                     posB=(pos_b, self.yt - 2 - i + self.scale * self.height), fc='#43CD80', ec='#43CD80',
                     arrowstyle='simple, head_width=5, head_length=5, tail_width=5', capstyle='butt')
                ax.add_artist(arrow)
            # 3.xorf text 
            if self.strand == '+':
                xorf_text_position = xorf[-1][1] * 1.01
            else:
                xorf_text_position = xorf[0][1] * 1.01
            plt.text(xorf_text_position, self.yt - 2 - i + self.scale * self.height, \
                     self.transcript_id + '.'+self.feature_type, fontsize=5)
        # arrow show the strand
        #if self.strand == '+':
        #    arrow = FancyArrowPatch(
        #        posA=(self.exon_location[0][0], 15), posB=(self.exon_location[-1][1], 15), fc=None, ec='#1f77b4',
        #        arrowstyle='simple, head_width=0.5, head_length=0.5, tail_width=0.2')
        #else:
        #    arrow = FancyArrowPatch(
        #        posA=(self.exon_location[-1][1], 15), posB=(self.exon_location[0][0], 15), fc=None, ec='#1f77b4',
        #        arrowstyle='simple, head_width=10, head_length=10, tail_width=5')

        ax.add_artist(arrow)
        # plt.axvline(x=self.exon_location[-1][1])
        plt.xlim(self.exon_location[0][0] - 100, self.exon_location[-1][1] + 500)
        plt.ylim(-2 * len(self.xorf_location), 10 + 2 * len(self.xorf_location))
        plt.axis('off')
        plt.savefig(f'{self.output}.pdf')
        plt.show()

#a = visual_transcript_genome('+', 'out', 'aaa',[[1,100],[200,300],[340,400],[440,550]], \
#                              [[50,100],[200,300],[340,400],[440,500]], \
#                              [[[20,100],[200,240]],[[20,100],[200,300],[340,360]],[[50,100],[200,300]]],'type')
#a.draw()
